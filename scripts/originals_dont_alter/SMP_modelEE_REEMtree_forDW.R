
#	This script models the SMP expert elicitation (EE) responses as functions of species traits and spatial context (threats, veg condition, etc).
# It then creates rasters of expected action benefits for each species x action combination.


library(mboost)
library(reshape2)
library(doParallel)
library(raster)
library(tcltk)

#number of parallel nodes (for doParrallel)
nnodes=2

# directory where EE data & spatial predictors stored
fdir_in <- "O:/SAN_Projects/SMP/SMP Sensitivity/InputData"

#where benefit rasters will be written - will have to change for each simulation
fdir_out<-"O:/SAN_Projects/SMP/SMP Sensitivity/BenefitRasters/"

load(file.path(fdir_in,"response.df.Rdf"))
load(file.path(fdir_in,"EEsitedata.Rdf"))

predictor.stack=stack(list.files(paste0(fdir_in,"/predictors"),pattern=".tif",full.names=T))

response.df$group=gsub("Bats","Mammals",response.df$group)

rsp0  = response.df


grpnms = c("Mammal","Bird","Amphibian","Reptile","Plant")
mstopvals=c(420,640,315,500,580)
mstopvals.na=c(200,300,50,100,100)

for (k in c(1:5)) {  #looping over taxon groups
  
  group = paste(grpnms[k],"s",sep="")
  
  responses0 = rsp0[rsp0$group==group,]
  
  load(file.path(fdir_in,paste0(grpnms[k],"Traits_EEModels.Rdf")))
  
  # assign the taxon specific trait data to 'traits'
  
  
  if (k==1) {
    traits = mamm.traits
      } else {
    if (k==2) {
      traits = bird.traits
    } else {
      if (k==3) {
        traits = frog.traits
      } else {
        if (k==4) {
          traits = reptile.traits
        } else {
          traits = plant.traits
        }
      }
    }
  }
  
  
  #Model Benefits...currently uses meanBen...i.e. experts best estimate
  
  responses = responses0[responses0$action != "Do Nothing",] 
  
  threats=sitedata[responses$ref.scenario,]
  
  #get scale, strata, climate change settings from scenario title
  
  scale=as.numeric(gsub(".*_","",gsub("ha_.*","",responses$title)))
  scen=gsub("_st.*","",responses$ref.scenario)
  strata=as.numeric(gsub("_sc.*","",gsub(".*_st","",responses$ref.scenario)))
  nocc=rep(0,length(strata))
  nocc[grep("_noCC_",responses$ref.scenario)]=1
  
  
  traitclass=unlist(lapply(traits,class))
  for(t in names(traitclass)[traitclass=="integer"]){traits[t]=as.numeric(traits[,t])}
  
  trait=traits[match(responses$TaxonCode,traits$TaxonCode),-c(1:3)] #get traits (removing ID variables)
  
  
  
  #model training data
  mdat=data.frame(y=as.numeric(responses$meanBen),action=responses$action,species=factor(responses$TaxonCode),expert=responses$expert,prob0=responses$refmean/10,tempscale=responses$tempscale,condition=responses$condition,scale,threats,strata,nocc,trait)
  
  mdat=mdat[is.finite(mdat$y),]
  
  mdat$species=factor(mdat$species,levels=c(levels(mdat$species),"newspp"))

  fmla=paste("y~btree(",paste(names(mdat)[!names(mdat)%in%c("y","expert")],collapse=","),")")
  fmla=as.formula(paste(fmla,"+ brandom(expert)"))  #+ brandom(scen) 
  
  modelBen<-mboost(fmla,data=mdat,control=boost_control(mstop=mstopvals[k],trace=T))
  
  
  ###Predict persistence probs conditional on No Action.
  
  responses = responses0[responses0$action == "Do Nothing",]
  
  threats=sitedata[responses$ref.scenario,]
  
  
  #get scale, strata, climate change settings from scenario title
  
  scale=as.numeric(gsub(".*_","",gsub("ha_.*","",responses$title)))
  scen=gsub("_st.*","",responses$ref.scenario)
  strata=as.numeric(gsub("_sc.*","",gsub(".*_st","",responses$ref.scenario)))
  nocc=rep(0,length(strata))
  nocc[grep("_noCC_",responses$ref.scenario)]=1
  
  trait=traits[match(responses$TaxonCode,traits$TaxonCode),-c(1:3)]
  
  #model training data
  
  mdat.na=data.frame(y=as.numeric(responses$refmean),action=responses$action,species=factor(responses$TaxonCode),expert=responses$expert,tempscale=responses$tempscale,condition=responses$condition,scale,threats,strata,nocc,trait)
  
  mdat.na=mdat.na[is.finite(mdat.na$y),]
  
  mdat.na$species=factor(mdat.na$species,levels=c(levels(mdat.na$species),"newspp"))
  #mdat.na$expert=factor(mdat.na$expert,levels=c(unique(mdat.na$expert),"AAA"))
  
  
  #build model
  fmla=paste("y~btree(",paste(names(mdat.na)[!names(mdat.na)%in%c("y","expert")],collapse=","),")")
  fmla.na=as.formula(paste(fmla,"+ brandom(expert)"))  #+ brandom(scen) 
  
  modelDN<-mboost(fmla.na,data=mdat.na,control=boost_control(mstop=mstopvals.na[k],trace=T))
  
  ##SAVE MODEL OBJECTS (if mapping predictions separately)
  
  #save(modelBen,file=file.path(fdir_out,"modelBen.Rmod"))
  #save(modelDN,file=file.path(fdir_out,"modelDN.Rmod"))
  
  
  ##MAP MODEL PREDICTIONS ACROSS EXTENT OF HDMS..... 

  
  savenames=c("Harvesting","Weeds","Rabbits","Dom.Grazing","Cats","Foxes","Fox.Cats","Graz.Weeds","Reveg","Deer","Fuelred","Goat","Horses","Miners","Phytoph","Pigs")
  actionnames=c("StopHarvesting","WeedControl","RabbitControl","StockExclusion","CatControl","FoxControl","FoxandCatControl","WeedControlandStockExclusion","Reveg","DeerControl","ReducedFRBfrequency","FeralGoatControl","FeralHorseControl","MinerControl","Phytophthoramanagement","FeralPigControl")
  
  vic=raster(file.path(fdir_in,"vicmask225vg.tif"))
  vicnv=raster(file.path(fdir_in,"nv225vg.tif"))
  
  threat.action=list(DeerControl="deer",FeralGoatControl="goat",FeralPigControl="pigs",
                     Phytophthoramanagement="phytophthora",RabbitControl="rabbit",ReducedFRBfrequency="fuelburning",
                     StockExclusion="grazing",StopHarvesting="harvesting",
                     WeedControl="transweed",WeedControlandStockExclusion=c("grazing","transweed"))
  
  traitnames=names(mdat)[names(mdat)%in%names(traits)]
  tscale=50
  scale=100
  nocc=0

  nvcells=which(values(vicnv)==1)
  
  predictor.vals=values(predictor.stack)
  
  if(group!="Plants"){
    hdmlist=list.files("O:/SAN_Projects/NPRINT4/HDMV52016/Tiff225/Masked/SMP_HDMS_Fauna100percentMaskExtant/",pattern=".tif$",full=T)
  }else{
    hdmlist=list.files("O:/SAN_Projects/NPRINT4/HDMV52016/Tiff225/Masked/SMP_HDMS/",pattern=".tif$",full=T)
    hdmlist2=list.files("O:/SAN_Projects/NPRINT4/HDMV52016/Tiff225/Masked/SMP_HDMS_MultiSppFlora",pattern=".tif$",full=T)
    hdmlist=c(hdmlist,hdmlist2)  
  }
  mapacts=levels(factor(responses0$action))
  mapacts=c("Do Nothing",mapacts[mapacts!="Do Nothing"])
  
  mapacts=mapacts[!mapacts%in%c("Phytophthora management","Reduced FRB frequency","Stop Harvesting")]
  
  
  #imput missing trait values (assuming model used can't handle NA's for prediction)
  
  if(!"TaxonGroup"%in%names(traits)){traits$TaxonGroup=1}
  
  traits.impt=traits  
  meanfun=function(x){
    if(is.numeric(x)){
      x=ifelse(is.na(x),mean(x,na.rm=T),x)
    }else{
      levs=levels(x)
      meanval=names(sort(table(x),decreasing=T))[1]
      x[is.na(x)]=meanval
    }
    return(x)
  }
  
  
  if(any(is.na(traits$TaxonGroup))){
    f2=factor(traits$TaxonGroup,levels=c(levels(traits$TaxonGroup),"MISSING"))
    f2[is.na(f2)]="MISSING"
    traits$TaxonGroup=f2
  }
  
  avgover=traits$TaxonGroup
  if(group=="Plants"){avgover=traits$florafauna_spp_lu_Family}
  for(c in 4:ncol(traits)){traits.impt[,c]=unlist(tapply(traits.impt[,c],avgover,meanfun))}
  
  NoRespSpp=traits$TaxonCode[traits$TaxonCode%in%responses$TaxonCode==F]
  
  doSpecies=traits$TaxonCode[1:5]
  
  cl <- makeCluster(nnodes)
  registerDoParallel(cl)
  foreach(sp=iter(doSpecies),.packages=c("raster","tcltk"),.errorhandling="stop")%dopar%{
    
    usehdm=grep(paste0("Spp",sp),hdmlist)
    if(length(usehdm)>0){
      hdm=raster(hdmlist[usehdm[1]])/100
      
      habvals=which(values(hdm)>0)
      
      pbsp <- tkProgressBar(paste("Mapping benefits for",sp), min=0, max=nlevels(mdat$action),label=paste(length(mapacts),"actions to do"))
     
       for(actn in mapacts){
        
        if(sp%in%NoRespSpp){possrange=quantile(mdat$y,c(0.1,.9))}else{
          possrange=range(mdat$y[mdat$species==sp &  mdat$action==actn],na.rm=T)/10
          if(sum(is.finite(possrange))<2){
            possrange=range(mdat$y[mdat$species==sp],na.rm=T)/10
          }
          if(sum(is.finite(possrange))<2){
            possrange=quantile(mdat$y,c(0.1,.9))
          }
        }
        
        usemod=modelBen
        usemdat=mdat
        
        if(actn=="Do Nothing"){
          usemod=modelDN
          usemdat=mdat.na
        }
        
        usevals=habvals
        aname=gsub(" ","",actn)
        
        if(aname%in%names(threat.action)){
          tvals=data.frame(predictor.vals[habvals,threat.action[[aname]]])
          usevals=habvals[which(apply(tvals,1,prod)>0)]
        }
        nvals=length(usevals)
        
        preds=vic*0
        #values(preds)=ifelse(values(vic)==1,0,NA)
        if(nvals>10){
          const=data.frame(species=factor(sp,levels=levels(usemdat$species)),action=factor(actn,levels=levels(usemdat$action)),tempscale=tscale,scale=scale,nocc=0)
          const=data.frame(const,traits.impt[match(sp,traits.impt$TaxonCode),traitnames])  
          
          if(sp%in%NoRespSpp){const$species=factor("newspp",levels=levels(usemdat$species))}
          
          newdat=data.frame(predictor.vals[usevals,],const)
          if(actn!="Do Nothing"){newdat$prob0=donothing[usevals]}
          
          
          ##THIS SPLITS RASTER UP INTO 200K cell chunks for processing...loops over chunks..PROB NOT EFFICIENT!
          nstarts=ceiling(nvals/200000)
          starts=seq(1,by=200000,length.out=nstarts)
          ends=c(starts[-1]-1,nvals)
          
          
          pb <- tkProgressBar(paste("Predicting",actn), min=0, max=nstarts,label=paste(nstarts,"subsets to do"))
          for(ind in 1:nstarts){
            newvals=predict(usemod,newdata=newdat[starts[ind]:ends[ind],],which=1)+usemod$offset
            values(preds)[usevals[starts[ind]:ends[ind]]]=newvals/10
            setTkProgressBar(pb, ind,label=paste(ind,"of",nstarts,"subsets done"))
            gc()
          }
          close(pb)
          if(actn=="Do Nothing"){
            values(preds)=ifelse(values(preds)<0,0,ifelse(values(preds)>1,1,values(preds)))
            donothing=values(preds)
            preds=preds*hdm
            DNpreds=values(preds)
          }
          
          
          
          if(actn!="Do Nothing"){
            preds=preds*hdm
            values(preds)=ifelse(values(preds)>(1-DNpreds),(1-DNpreds),ifelse(values(preds)<(-DNpreds),-DNpreds,values(preds)))
          }
          
        }
        
        shortname=savenames[match(aname,actionnames)]
        
        tifname=paste0("Spp",sp,"delta_",shortname,".tif")
        if(actn=="Do Nothing"){tifname=paste0("Spp",sp,"_DoNothing.tif")}
        
        writeRaster(preds,file=file.path(fdir_out,tifname),format="GTiff",overwrite=T)
        
        setTkProgressBar(pbsp, match(actn,mapacts),label=paste(actn, "done"))  
        
      }
      close(pbsp)
      gc()
    }
    if(!exists("pb0")) pb0 <- tkProgressBar("Mapping Species Benefits", min=0, max=length(doSpecies))
    setTkProgressBar(pb0, match(sp,doSpecies))
  }
  stopCluster(cl)
  gc()
  
  
}















