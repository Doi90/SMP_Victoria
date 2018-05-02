

## Create 'delta' staks for each species, and calculated weighted sum of delta values (benefits)
## This is the basis of stage 1 prioritization. Nspecies X Nactions rasters built- long processing times!!

collateBenefits=function(do.species,# list of taxoncodes of the species to include
                         sppInfo, # table with species weights, dist. sums etc
                         BenefitDir, # directory where benefit rasters are stored
                         rasterStack.folder, #folder where (smoothed) benefit raster stacks will be stored (will be created if doesn't exist)
                         actions=actions,#action names
                         nvactions=actions[-grep("Reveg",actions)],# actions that apply to native veg (generally all but reveg)
                         AddSpecies=F,addTo=NULL,  #TRUE if adding species to existing totben.vals table (addTo = full path to that R object)
                         nnodes=12, #number of parallel nodes to use...
                         mask=raster("O:/SAN_Projects/SMP/vicmask225.tif"), # analysis mask (cells with value >0 included in analysis)
                         vegmask=raster("O:/SAN_Projects/SMP/vicNV225_binary.tif"), # native veg mask (> 0 for native veg)
                         costs, 
                         smoothbenefits=F, # smooth benefits over minimum area?
                         mincells=NULL,# names vector of min. no. of pixels over which an action should be performed.
                         onParError="stop", #"stop" for debugging (will cause break if something missing)
                         rastempdir="E:TEMP/",
                         maxgpsize=400,deletemp=T
){# function calcs start here  

  rasterOptions(tmpdir=rastempdir)
  
  BF=function(R,x=0.25,y=1,T=1,w1=1,w2=1){ifelse(R<=T,w1*((R/T)^x),w1+w2*(((R-T)/(1-T))^y))}
  na0=function(x){ifelse(is.na(x),0,x)}
  
  mask=mask>0
  clipvals=ifelse(is.na(values(mask)),0,1)
  viccells=which(clipvals==1)
  outcells=which(clipvals!=1)
  nvcells=which(values(vegmask>0)==1) 
  
  
  setwd(rasterStack.folder)
  totbens=mask-1
  for(a in 2:length(actions)){totbens=addLayer(totbens,(mask-1))}
  names(totbens)=actions
  totben.zero=values(totbens)[viccells,]
  rm(totbens)
  
  if(AddSpecies){
    load(addTo)
    load(gsub("totben.vals","weightings",addTo))
    #isdone=gsub("_smdeltas.Rstack","",list.files(pattern="_smdeltas.Rstack"))
    #do.species=species[species%in%isdone==FALSE]
  }
  if(AddSpecies==FALSE){totben.vals=totben.zero; weightings=NULL}
  
  gc()
  
  #set up batches of  maxgpsize spp to loop thru...to avoid memory leek problem...(find a better fix!)
  if(length(do.species)< maxgpsize){ngps=1;gpid=rep(1,length(do.species))}else{
    ngps=ceiling(length(do.species)/ maxgpsize)
    gpid=rep(1, maxgpsize)
    if(ngps>2){
    for(gpn in 2:(ngps-1)){
      gpid=c(gpid,rep(gpn, maxgpsize))
    }}
    gpid=c(gpid,rep(ngps,length(do.species)-length(gpid)))
  }
  
  anyna=function(x){any(is.na(x))}
  
  combineFun=function(list1,list2){
    list(list1[[1]]+list2[[1]],rbind(list1[[2]],list2[[2]]))
  }
  
  #BenefitDir=settings["BenefitRastersDirectory","Value"]
  benlayers=list.files(BenefitDir,pattern=".tif$",full=T)
  dudly=grep("delta_DoNothing",benlayers)
  if(length(dudly)>0){
    benlayers=benlayers[-grep("delta_DoNothing",benlayers)]
  }
  
  rm(pb0)
  
  #nvactions=actions[-grep("Reveg",actions)]
  for(gp in 1:ngps){
    
    cl <- makeCluster(nnodes)
    registerDoParallel(cl)
    start.time=Sys.time()

    smoothit<-foreach(sp=iter(do.species[which(gpid==gp)]),.packages=c("raster","tcltk"),.combine=combineFun,.errorhandling=onParError)%dopar%{
      spname=rownames(sppInfo)[match(sp,sppInfo$taxoncode)]
      nameroot=paste0("Spp",sp,"delta_")
      spbens=benlayers[grep(nameroot,benlayers)]
      smdeltas.nas0=totben.zero
      maxben_sdmax=data.frame(maxben=NA,sum.delta.max=NA,minBF=NA,minprop=NA,maxprop=NA,sumreveg=NA)
      rownames(maxben_sdmax)=sp
      if(length(spbens)>0){
        deltas=stack(spbens,quick=T)
        #values(deltas)[outcells,]=NA
        names(deltas)=gsub(nameroot,"",names(deltas))

        smdeltas=deltas
        rm(deltas)
        for(actn in nvactions){
          if(actn%in%names(smdeltas)){
            values(smdeltas[[actn]])[viccells]=ifelse(is.na(values(costs[[actn]])[viccells]),0,values(smdeltas[[actn]])[viccells])
          }}
   
        spdonot=benlayers[grep(paste0("Spp",sp,"_DoNothing.tif"),benlayers)]
        donothing=raster(spdonot)
        dnvals=values(donothing)
        sumdonot=sum(dnvals[viccells],na.rm=T)
        maxdelta.vals=values(max(smdeltas[[nvactions]],na.rm=T)) #reveg should be included in this calc IF confident about reveg benefits
        maxhab=dnvals+maxdelta.vals

        minprop=sumdonot/sppInfo[spname,"sum1750cells"]
        maxprop=sum(maxhab,na.rm=T)/sppInfo[spname,"sum1750cells"]
        minBF=BF(min(minprop,1),x=sppInfo[spname,"Vx"])
        maxben=BF(min(maxprop,1),x=sppInfo[spname,"Vx"])-minBF
        if(!is.finite(maxben)){maxben=1}

        sum.delta.max=sum(maxdelta.vals,na.rm=T)
        if(sum.delta.max==0){sum.delta.max=1}
        sumreveg=0
        if("Reveg"%in%names(smdeltas)){sumreveg=sum(values(smdeltas[["Reveg"]]),na.rm=T)}
        maxben_sdmax=data.frame(maxben,sum.delta.max,minBF,minprop,maxprop,sumreveg)
        rownames(maxben_sdmax)=sp
      #
        if(smoothbenefits){
          for(actn in actions){
            if(mincells[actn]>1){
              nbv=round(sqrt(mincells[actn]))
              if(nbv%%2==0){nbv=nbv+1}
              smdeltas[[actn]]=focal(smdeltas[[actn]], w=matrix(1,nbv,nbv), fun=mean,na.rm=T,pad=F)
              values(smdeltas[[actn]])=ifelse(is.na(values(costs[[actn]])),NA,1)*values(smdeltas[[actn]])
              gc()
            }
          }}
      #
        save(smdeltas,file=file.path(rasterStack.folder,paste0("Spp",sp,"_smdeltas.Rstack")))
        smdeltas.nas0=totben.zero
        if(sum(apply(values(smdeltas)[viccells,],2,anyna))>0){
          values(smdeltas)[viccells,]=apply(values(smdeltas)[viccells,],2,na0)
        }
        if(is.finite(sppInfo[spname,"weight"]*maxben/sum.delta.max)){
          smdeltas.nas0[,names(smdeltas)]=values(smdeltas)[viccells,]*sppInfo[spname,"weight"]*maxben/sum.delta.max
        }
        rm(smdeltas)
        gc()
        if(!exists("pb0")) pb0 <- tkProgressBar("Building benefit outputs", min=0, max=length(do.species))
        setTkProgressBar(pb0, match(sp,do.species))
      }
      list(smdeltas.nas0,maxben_sdmax)
    }
    
    end.time=Sys.time()
    stopCluster(cl)
    
    weightings=rbind(weightings,smoothit[[2]])
    totben.vals=totben.vals+smoothit[[1]]
    #gotem=rownames(smoothit[[2]])[rownames(smoothit[[2]])%in%sppInfo$taxoncode]
    #sppInfo[match(as.numeric(gotem),sppInfo$taxoncode),names(smoothit[[2]])]=smoothit[[2]][gotem,]
    #save(smoothit[[2]],file=paste0("smoothit_",gp,".R"))
    save(totben.vals,file=paste0("totben.vals_",gp,".Rmat"))
    save(weightings,file=paste0("weightings_",gp,".Rmat"))
    gc()
  }
  
  tbenfile.name=file.path(rasterStack.folder,"totben.vals.Rmat")
  weightfile.name=file.path(rasterStack.folder,"weightings.Rmat")
  save(totben.vals,file=tbenfile.name)
  save(weightings,file=weightfile.name)
  
  message(cat("Done compiling summed benefits.\n",paste("Matrix of weighted benefit sums saved at:",tbenfile.name,sep="\n")))
  message(cat("\n",paste("Weighting factors for Zonation saved at:",weightfile.name,sep="\n")))
  

  if(deletemp){
  tempfiles=c(list.files(pattern="totben.vals_",full.names = T),list.files(pattern="weightings_",full.names=T))
  file.remove(tempfiles)
  }
  
  return(list(totben.file=tbenfile.name,zweights.file=weightfile.name))
}


getBestActions=function(benefits,costs,clearing.probs,nogo.combos,clearing.ind=1,all=999,outpath=output.folder,mask=raster("O:/SAN_Projects/SMP/vicmask225.tif"),viccells=NULL,scen.name="scen",is.shiny=F)
{
  if(is.null(viccells)){
    mask=mask>0
    clipvals=ifelse(is.na(values(mask)),0,1)
    viccells=which(clipvals==1)
  }
  for(r in 1:nrow(nogo.combos)){nogo.combos[r,]=match(nogo.combos[r,],colnames(benefits))}
  nogo.combos=apply(nogo.combos,2,as.numeric)  
  
  if(is.shiny){
  withProgress(min=0,max=nrow(benefits),value=0,message="Finding cost efficient actions", {
  maxs=sapply(1:nrow(benefits), function(i) argmax.comb.shiny(progind=i,ben=benefits[i,], cost=costs[i,],clearing.prob=clearing.probs[i],nogos=nogo.combos,clearing.ind=clearing.ind,all=all))
  })
  }else{
    maxs=sapply(1:nrow(benefits), function(i) argmax.comb(ben=benefits[i,], cost=costs[i,],clearing.prob=clearing.probs[i],nogos=nogo.combos,clearing.ind=clearing.ind,all=all))
  }
  
  maxs=data.frame(t(maxs))
  save(maxs,file=file.path(outpath,paste0("maxcombos_",scen.name,".Rmat")))
  actcodes=ifelse(unlist(maxs$bestcol)==999,ncol(benefits)+as.numeric(as.factor(unlist(maxs$bestid))),unlist(maxs$bestcol))
  bestact=mask
  values(bestact)[viccells]=actcodes
  
  lookup=unique(data.frame(actcodes,combo=unlist(maxs$bestid)))
  lookup=lookup[order(lookup$actcodes),]
  lookup=data.frame(lookup,t(data.frame(lapply(as.vector(lookup$combo),namethem,actnames=colnames(benefits)))))
  lookup$actions=unlist(lapply(1:nrow(lookup),function(i) paste(na.omit(unlist(as.vector(lookup[i,3:5]))),collapse="+")))
  lookup=lookup[lookup$actions!="",]  
  
  LUname=file.path(outpath,paste0(scen.name,"_BAlookup.csv"))
  write.csv(lookup,LUname)
  
  bestact=bestact*mask
  bestact=ratify(bestact)
  levs=levels(bestact)[[1]]
  
  #save bestaction rasters in single and multi-layer formats
  baRname=file.path(outpath,paste0("bestaction_",scen.name,".Rras"))
  baTifName=file.path(outpath,paste0("bestaction_",scen.name,".tif"))
  baCostTifName=file.path(outpath,paste0("bestaction_",scen.name,"_costs.tif"))
  
  save(bestact,file=baRname)
  writeRaster(bestact,file=baTifName,format="GTiff",overwrite=T)
  
  ba.stack=addLayer(mask,addLayer(mask,mask))
  values(ba.stack[[1]])=match(lookup$X1[match(values(bestact),lookup$actcodes)],colnames(benefits))
  values(ba.stack[[2]])=match(lookup$X2[match(values(bestact),lookup$actcodes)],colnames(benefits))
  values(ba.stack[[3]])=match(lookup$X3[match(values(bestact),lookup$actcodes)],colnames(benefits))
  names(ba.stack)=paste0("bestact",1:3)
  baTif3name=file.path(outpath,paste0("SMPbestaction_3band_",scen.name,".tif"))
  baTif3.lu.name=file.path(outpath,paste0("SMPbestaction_3band_",scen.name,"_lookup.csv"))
  writeRaster(ba.stack,file=baTif3name,overwrite=T)
  
  balookup2=data.frame(id=c(1:ncol(benefits)),action=colnames(benefits))
  write.csv(balookup2,file=baTif3.lu.name)
  
  costs.best=mask
  values(costs.best)[viccells]=as.numeric(maxs$cost)
  writeRaster(costs.best,file=baCostTifName,format="GTiff",overwrite=T)
  
  
  message(cat("Done getting best local actions.\n","Result rasters saved to:"))
  message(cat(baRname,"(R raster)","\n",baTifName,"(GeoTiff)","\n",baTif3name,"(multi-band GeoTiff)","\n",LUname,"(lookup table)","\n",baCostTifName,"(GeoTiff)"))
  
  return(list(raster=baRname,lookup=lookup,cost.raster=baCostTifName))
  
}



getBAbenefits=function(bestact,lookup,rasterStack.folder,features,species,BenefitDir,actions,
                       nvcells,species.done=NULL,AddSpecies=FALSE,Zwd,mask,nnodes=12,outpath,onParError="remove"){
  
  stackSelect.fast=function(Stack,index,base=clip,nato0=TRUE){
    outras=base*0
    stackvals=values(Stack)
    for(ind in unique(na.omit(index))){
      newvals=stackvals[which(index==ind),ind]
      if(nato0){newvals=ifelse(is.na(newvals),0,newvals)}
      values(outras)[which(index==ind)]=newvals
    }
    return(outras)
  }
  
  noact.sums=features$sum1750cells*features$minprop
  weights.adj=nzcount=smdelta.best.sums=NULL
  
  new.species=species
  
  if(AddSpecies){
    load(paste0(outpath,"/base.sums.Rvect"))
    load(paste0(outpath,"/noact.sums.Rvect"))
    load(paste0(outpath,"/weights.adj.Rvect"))
    load(paste0(outpath,"/nzcount.Rvect"))
    load(paste0(outpath,"/base.std.Rstack"))
    done.species=species[species%in%species.done]
    new.species=species[species%in%species.done==FALSE]
  }
  
  new.slots=rep(NA,length(new.species))
  names(new.slots)=new.species
  weights.adj=c(weights.adj,new.slots)
  nzcount=c(nzcount,new.slots)
  
  setwd(rasterStack.folder)
  
  None=mask-1
  
  best.index=values(bestact)
  best.index.comb=matrix(nrow=length(best.index),ncol=3)
  best.index.comb[,1]=match(lookup$X1,actions)[match(best.index,lookup$actcodes)]
  best.index.comb[,2]=match(lookup$X2,actions)[match(best.index,lookup$actcodes)]
  best.index.comb[,3]=match(lookup$X3,actions)[match(best.index,lookup$actcodes)]
  
  best.index.comb0=best.index.comb

  #balookup2=data.frame(id=c(1:length(actions)),action=actions)
  #write.csv(balookup2,file=file.path(outpath,"BAlookup_multiband.csv"))
  
  countnz=function(x){sum(ifelse(na.omit(x)>0,1,0))}
  
  base.sums=features$sum1750cells
  if(is.null(base.sums)){base.sums=features$sum1750ha}
  names(base.sums)=features$taxoncode
  
  nnodes=nnodes
  cl <- makeCluster(nnodes,outfiles="")
  registerDoParallel(cl)
  start.time=Sys.time()
  
  getbestvals<-foreach(sp=iter(species),.packages=c("raster","tcltk"),.combine=rbind,.errorhandling=onParError)%dopar%{
    
    load(paste0("Spp",sp,"_smdeltas.Rstack"))
    donot=raster(file.path(BenefitDir,paste0("Spp",sp,"_DoNothing.tif")))
    #values(donot)[viccells]=ifelse(is.na(values(donot)[viccells]),0,values(donot)[viccells])
    missing.bens=actions[actions%in%names(smdeltas)==F]
    if(length(missing.bens)>0){
      for(mb in missing.bens){
        smdeltas=addLayer(smdeltas,None)
        names(smdeltas)[nlayers(smdeltas)]=mb
      }
    }
    smdeltas=smdeltas[[actions]]
    
    if(sp%in%new.species){
      weights.adj.sp=round(features[match(sp,features$taxoncode),"Zweight"],3) 
    }
 
    save(smdeltas,file="smdeltas.Rstack")
    best.smdelta1=stackSelect.fast(smdeltas,best.index.comb[,1],base=mask)
    best.smdelta2=stackSelect.fast(smdeltas,best.index.comb[,2],base=mask)
    best.smdelta3=stackSelect.fast(smdeltas,best.index.comb[,3],base=mask)
    
     best.smdelta=best.smdelta1+best.smdelta2+best.smdelta3
    fincells=which(is.finite(values(best.smdelta)))
    values(best.smdelta)[fincells]=ifelse(values(best.smdelta)[fincells]>1-values(donot)[fincells],1-values(donot)[fincells],values(best.smdelta)[fincells])
    #values(best.smdelta)[viccells]=ifelse(is.na(values(best.smdelta)[viccells]),0,values(best.smdelta)[viccells])
     sumneg=NA
    mind=minValue(best.smdelta)
    if(mind<0){
      negvals=which(values(best.smdelta)<0)
      absnegvals=abs(values(best.smdelta)[negvals])
      values(best.smdelta)[negvals]=0
      if(mind<(-0.01) & length(negvals)>10){
        negras=mask-1
        values(negras)[negvals]=absnegvals
        writeRaster(negras,file=paste0(Zwd,"/inputs/Spp",sp,"_negdelta.tif"),format="GTiff",overwrite=T)
        sumneg=sum(absnegvals,na.rm=T)
      }}
    
    writeRaster(best.smdelta,file=paste0(Zwd,"/inputs/Spp",sp,"_delta.tif"),format="GTiff",overwrite=T)
    sumpos=sum(values(best.smdelta),na.rm=T)
    
    smdelta.best.sums.sp=sum(values(best.smdelta),na.rm=T)
    
 
    weight.row=data.frame(weights.adj.sp,smdelta.best.sums.sp,sumpos,sumneg)
    rownames(weight.row)=sp
    
    if(!exists("pb0")) pb0 <- tkProgressBar("Making benefit layers for Zonation", min=0, max=length(species))
    setTkProgressBar(pb0, match(sp,species))
    gc()
    weight.row
  }
  
  end.time=Sys.time()
  stopCluster(cl)
  
  weight.rows=getbestvals
  
  # noact.sums[rownames(weight.rows)]=weight.rows[,"noact.sums.sp"]
  # weights.adj[rownames(weight.rows)]=weight.rows[,"weights.adj.sp"]
  # nzcount[rownames(weight.rows)]=weight.rows[,"nzcount.sp"]
  # smdelta.best.sums[rownames(weight.rows)]=weight.rows[,"smdelta.best.sums.sp"]
  sumnegras=weight.rows$sumneg
  sumposras=weight.rows$sumpos
  names(sumnegras)=names(sumposras)=rownames(weight.rows)
  
  weight.rows$qweight=weight.rows$smdelta.best.sums.sp/base.sums[rownames(weight.rows)]
  
  missing=names(which(is.na(weight.rows$weights.adj.sp)))
  
  weightMatName=file.path(outpath,"weight.rows.R")
  
  save(weight.rows,file=weightMatName)
  
  
  message(cat("Done creating benefit rasters for Zonation.\n","Rasters saved to:",paste0(Zwd,"/inputs/")))
  message(cat("weight rows saved to:",weightMatName))
  
  return(weight.rows)
  
  
}



makeZfiles=function(Zwd,features,weight.rows,species=species,
                  removal.rule=2,
                  warp=10000,
                  edge.removal=1,
                  edge.points=10000,
                  scen.name="scen",
                  cost.best.path=NULL,
                  firecost=NA,
                  zigdir="E:/SMPv1.0/Zonation"
                  
){

setwd(Zwd)

if(!is.null(cost.best.path)){file.copy(cost.best.path,file.path(Zwd,"/inputs/SMP_costbest.tif"))}

datname=file.path(Zwd,paste0("SMP_zonation_settings_",scen.name,".dat"))

sink(file=datname)

cat("[Settings]","\n")
cat("removal rule =",removal.rule,"\n")
cat("warp factor =",warp,"\n")
cat("edge removal =",edge.removal,"\n")
cat("add edge points =",edge.points,"\n","\n")

cat("use cost = 1","\n")
cat("cost file =",paste0(Zwd,"/inputs/SMP_costbest.tif"),"\n")

firecost=NA
if(!is.na(firecost)){file.copy(firecost,to=paste0(Zwd,"/inputs/firecost.tif"))}

sink()


#Write tifs and Zonation features list file
Zdeltafiles=list.files(paste0(Zwd,"/Inputs/"),pattern="_delta.tif$",full=T)
spplistname=file.path(Zwd,paste0("SMP_species_list_",scen.name,".spp"))
sink(file=spplistname)

negwgt.adj=-weight.rows$qweight*weight.rows$sumneg/weight.rows$sumpos 
names(negwgt.adj)=rownames(weight.rows)

for(sp in species){
  filename=Zdeltafiles[grep(sp,Zdeltafiles)]
  if(is.finite(features[match(sp,features$taxoncode),"Zweight"]) & length(filename)==1){
    cat(weight.rows[toString(sp),"qweight"]*features[match(sp,features$taxoncode),"Zweight"],0.001,1,1,features[match(sp,features$taxoncode),"expon"],filename,"\n")  #removed from weight '*smdelta.best.sums[sp]'
    if(file.exists(gsub("_delta","_negdelta",filename))){
      cat(negwgt.adj[toString(sp)]*features[match(sp,features$taxoncode),"Zweight"],0.001,1,1,features[match(sp,features$taxoncode),"expon"],gsub("_delta","_negdelta",filename),"\n") 
    }
  }
}

if(!is.na(firecost)){
  cat(-1,0.001,1,1,1,paste0(Zwd,"/Inputs/firecost.tif"),"\n")
}
sink()

#Write Zonation .bat file
outname=file.path(Zwd,paste0("SMP_",scen.name,".txt"))
Zcall=paste("call",file.path(zigdir,"zig4.exe -r"),datname ,spplistname, outname, 0, 0, 1.0, 0)

sink(file=paste0("SMP_runZ_",scen.name,".bat"))
cat(Zcall,"\n")
sink()

message(cat("Done writing run files for Zonation.\n","Files saved in:",Zwd))


return(list(batfile=file.path(Zwd,paste0("SMP_runZ_",scen.name,".bat")),datfile=file.path(Zwd,datname),
            spplist=file.path(Zwd,spplistname),outname=outname))

}


# files to re-run ZOnation with limited reveg (reveg.lim %)
Zreveglim=function(Zwd,Zfiles=Zfiles,reveg.lim=2,zigdir="E:/SMPv1.0/Zonation",
                  revegcands=raster("O:/SAN_Projects/NPRINT4/Masks/Not_NV_Not_Reservoi_or_Urban.tif"),
                  nv=raster("O:/SAN_Projects/SMP/vicNV225_binary.tif")
){
setwd(Zwd)
rankfile=list.files(pattern="rank.compressed")
scen.name=gsub(".*\\/","",gsub(".txt","",Zfiles$outname))
rankfile=rankfile[grep(scen.name,rankfile)][1]

zrank=raster(rankfile)

#create analysis masks that excludes 90% & 95% of reveg 

rvranks=zrank*revegcands
rankrv=1-values(rvranks)
topXthresh=quantile(rankrv,reveg.lim/100,na.rm=T)

zrvgmaskX=nv
values(zrvgmaskX)=ifelse(is.na(values(nv)),ifelse(rankrv<=topXthresh,1,NA),values(nv))

zmasknameX=paste0(Zwd,"/inputs/zmask",reveg.lim,".tif")

writeRaster(zrvgmaskX,file=zmasknameX,format="GTiff",overwrite=T)



#append .dat file

newdat=gsub(".dat",paste0("_reveglim",reveg.lim,".dat"),Zfiles$datfile)
file.copy(Zfiles$datfile,newdat)

sink(newdat,append=T)
cat("area mask file = ",zmasknameX,"\n")
sink()

#append .bat file

newbat=gsub(".bat",paste0("_reveglim",reveg.lim,".bat"),Zfiles$batfile)

newout=gsub(".txt",paste0("_reveglim",reveg.lim,".txt"),Zfiles$outname)

ZcallX=paste("call",file.path(zigdir,"zig4.exe -r"),newdat ,Zfiles$spplist, newout, 0, 0, 1.0, 0)

sink(file=newbat)
cat(ZcallX,"\n")
sink()

Zfiles$outname=newout
Zfiles$datfile=newdat
Zfiles$batfile=newbat

message(cat("Done writing reveg. limited run files for Zonation.\n","Files saved in:",Zwd))


return(Zfiles)
}


#######NOW RUN IN ZONATION########
####### MAIN ANALYSIS STOPS HERE ######### 






#Make polygons for display best action prioprity...

polyfy=function(bestact,zrank,lookup.path,nzgps=5,scen.name,outpath=getwd(),replace=T,dtype="INT2U"){
  
lookup=read.csv(lookup.path)
  
predators=c("Foxes","Cats","Fox.Cats")
herbivores=c("Rabbits","Deer","Goat","Dom.Grazing","Graz.Weeds","Horses")
weeds=c("Weeds","Graz.Weeds")
other=c("Pigs","Miners")

lookup.predators=lookup$X1%in%predators | lookup$X2%in%predators | lookup$X3%in%predators
lookup.herbivors=lookup$X1%in%herbivores | lookup$X2%in%herbivores | lookup$X3%in%herbivores

lookup.clearing=lookup$X1=="Clearing" | lookup$X2=="Clearing" | lookup$X3=="Clearing"
lookup.clearing=ifelse(is.na(lookup.clearing),FALSE,lookup.clearing)
lookup.weeds=lookup$X1%in%weeds | lookup$X2%in%weeds | lookup$X3%in%weeds

lookup.other=lookup$X1%in%other | lookup$X2%in%other | lookup$X3%in%other

lookup.reveg=lookup$X1=="Reveg" | lookup$X2=="Reveg" | lookup$X3=="Reveg"
lookup.reveg=ifelse(is.na(lookup.reveg),FALSE,lookup.reveg)

lookup$action.class=ifelse(lookup.reveg,"Reveg",ifelse(lookup.clearing,"Clearing",
                    ifelse(lookup.predators & !lookup.herbivors & !lookup.weeds,"predators",
                    ifelse(lookup.herbivors & !lookup.predators & !lookup.weeds,"herbivores",
                    ifelse(lookup.weeds & !lookup.predators & !lookup.herbivors,"weeds",
                    ifelse(lookup.herbivors & lookup.predators & !lookup.weeds,"preds.herbvrs",
                    ifelse(lookup.herbivors & !lookup.predators & lookup.weeds,"herbvrs.weeds",
                    ifelse(lookup.herbivors & lookup.predators & lookup.weeds,"preds.herbvrs.weeds",
                    ifelse(lookup.predators & !lookup.herbivors & lookup.weeds,"preds.weeds",
                    ifelse(lookup.other,"other",NA))))))))))

write.csv(lookup,lookup.path)

bestact.poly=rasterToPolygons(bestact)  
  
zband=zrank
values(zband)=ceiling(values(zrank)*10)
writeRaster(zband,file=file.path(outpath,paste0(scen.name,"_zrank_deciles.tif")),format="GTiff",overwrite=replace,datatype=dtype)  
zband5=zband
values(zband5)=ceiling(values(zband)/2)
writeRaster(zband5,file=file.path(outpath,paste0(scen.name,"_zrank_quintiles.tif")),format="GTiff",overwrite=replace,datatype=dtype)  
 
writeRaster(bestact,file=file.path(outpath,paste0(scen.name,"_bestactions.tif")),format="GTiff",overwrite=replace,datatype=dtype) 

##converting multiband to 3 separate tifs for ARCMAP...

for(l in 1:3){
  writeRaster(multi[[l]],file=file.path(outpath,paste0(scen.name,"_bestaction",c("A","B","C")[l],".tif")),format="GTiff",overwrite=replace,datatype=dtype) 
}
}


#smooth best actions...not currently used
#SmoothBestAction=settings["SmoothBestAction","Value"]

smoothBestActions=function(){
  
  if(SmoothBestAction){
    add.action=NULL
    alltrue=function(x){if(length(which(x==TRUE))==length(x)){TRUE}else{FALSE}}
    
    for(actn in levs$action[levs$action!="None"]){
      actn.ras=bestact
      actn.no=levs$ID[match(actn,levs$action)]
      values(actn.ras)=ifelse(values(bestact)==actn.no,1,0)
      nbv=max(5,sqrt(mincells[actn]))
      if(nbv%%2==0){nbv=nbv+1}
      actn.mean=focal(actn.ras, w=matrix(1,nbv,nbv), fun=mean,na.rm=T,pad=T)
      #values(actn.mean)=values(actn.mean)*ifelse(is.na(values(costs[[actn]])),NA,1)
      addit=ifelse(values(actn.mean)>=0.5 & values(actn.ras)==0,values(actn.mean),0)
      add.action=cbind(add.action,addit)
      colnames(add.action)[ncol(add.action)]=actn
      rm(actn.mean)
    }
    changeto=as.list(c("None",actions)[maxs+1])
    anychange=which(apply(add.action,1,sum,na.rm=T)>0)
    for(changeno in anychange){
      orig=changeto[[changeno]]
      orig.threats=threats[which(action.info[orig,threats]==1)]
      if(length(orig.threats)==1){possnew=actions[which(action.info[,orig.threats]==1)]}
      if(length(orig.threats)>1){possnew=names(which(apply(action.info[,orig.threats]==1,1,alltrue)))}
      possnew=intersect(names(which(add.action[changeno,]>0)),possnew)
      if(length(possnew)>1){possnew=names(which(add.action[changeno,possnew]==max(add.action[changeno,possnew])))}
      if(length(possnew)==0){possnew=names(which(add.action[changeno,]==max(add.action[changeno,])))}
      changeto[[changeno]]=possnew
    }
    
    values(bestact)[anychange]=match(unlist(changeto)[anychange],actions)+1
    
    rm(changeto);rm(anychange)
    
    plot(bestact,rainbow(length(actions)))
  }
}

#end smoothing




viewresults<-function(){

####Species performance curves
library(reshape2)
library(ggplot2)

curves.files=list.files(pattern=".curves.txt")
curves.file=curves.files[grep(date,curves.files)]

curves=read.table(curves.file)#,header=T)

##species info.
info.files=list.files(pattern=".features_info.txt")
info.file=info.files[grep(date,curves.file)]
feature.info=read.table(info.file,skip=1,header=T)
mapfiles.to.use=feature.info$MapFileName[feature.info$Weight>0 & feature.info$distribution.sum>0]

Zspecies=gsub(".tif","",gsub(paste0(Zwd,"/Inputs/"),"",mapfiles.to.use))
Zspecies=Zspecies[Zspecies%in%species]

names(curves)=c("Prop_landscape_lost","Cost","min_prop_rem","ave_prop_rem","W_prop_rem","ext1","ext2",gsub(".tif","",gsub(paste0(Zwd,"/Inputs/"),"",feature.info$MapFileName)))

curves.habvals=data.frame(t(t(curves[Zspecies])*smdelta.best.sums[Zspecies]+noact.sums[Zspecies]))
curves.deltas=data.frame(t(t(curves[Zspecies])*smdelta.best.sums[Zspecies]))
curves.deltas.relnow=sweep(curves.habvals,2,current.sums[Zspecies],"-")
curves.Rvals=sweep(curves.habvals,2,base.sums[Zspecies],"/")
curves.deltas.prop=sweep(curves.deltas,2,base.sums[Zspecies],"/")

curves.Rvals.relnow=sweep(curves.habvals,2,current.sums[Zspecies],"/")
curves.deltas.prop.relnow=sweep(curves.deltas.relnow,2,current.sums[Zspecies],"/")

curves.habvals=data.frame(apply(curves.habvals,2,function(x) {ifelse(x<0,0,x)}))
curves.Rvals=data.frame(apply(curves.Rvals,2,function(x) {ifelse(x<0,0,ifelse(x>1,1,x))}))
curves.deltas.prop=data.frame(apply(curves.deltas.prop,2,function(x) {ifelse(x<0,0,ifelse(x>1,1,x))}))

curves.Rvals.relnow=data.frame(apply(curves.Rvals.relnow,2,function(x) {ifelse(x<0,0,ifelse(x>5,5,x))}))
curves.deltas.prop.relnow=data.frame(apply(curves.deltas.prop.relnow,2,function(x) {ifelse(x<(-5),-5,ifelse(x>5,5,x))}))



perfplot.habvals=data.frame(curves$Cost,melt(curves.habvals))
names(perfplot.habvals)=c("Cost","Species","Habitat")

pp.habvals=ggplot(data=perfplot.habvals)+geom_line(aes(y=Habitat,x=Cost,group=Species,col=Species))
pp.habvals+ggtitle("Total habitat")

perfplot.deltas=data.frame(curves$Cost,melt(curves.deltas))
names(perfplot.deltas)=c("Cost","Species","Benefit")

pp.habvals=ggplot(data=perfplot.deltas)+geom_line(aes(y=Benefit,x=Cost,group=Species,col=Species))
pp.habvals+ggtitle("Benefits")#+ylim(0,)

perfplot.Rvals=data.frame(curves$Cost,melt(curves.Rvals))
names(perfplot.Rvals)=c("Cost","Species","proportion")

pp.Rvals=ggplot(data=perfplot.Rvals)+geom_line(aes(y=proportion,x=Cost,group=Species,col=Species))
pp.Rvals+ggtitle("Proportional Habitat")

curves.Vvals=curves.Rvals
for(sp in Zspecies){curves.Vvals[sp]=BF(curves.Rvals[sp],x=features[sp,"Vx"],T=features[sp,"VT"])}

perfplot.Rquants=data.frame(curves$Cost,(1-curves[,1])*totarea,t(apply(curves.Vvals,1,quantile,probs=c(0,0.5,0.1,0.25,0.5,0.75,0.9,0.95,1))))
names(perfplot.Rquants)=c("Cost","Area",c(0,5,10,25,50,75,90,95,100))
pplot.rquants=melt(perfplot.Rquants,id.vars=c("Cost","Area"))
names(pplot.rquants)=gsub("variable","Quantile",names(pplot.rquants))
pplot.rquants$Cost=pplot.rquants$Cost/1000000
pp.Rquant=ggplot(data=pplot.rquants)+geom_line(aes(y=value,x=Cost,group=Quantile,col=Quantile),lwd=2)
pp.Rquant+theme_bw(base_size=16)+xlab("Annual investment ($ Million)")+ylab("Persistence Probability")+xlim(c(0,20))

pp.Rquant+theme_bw(base_size=16)+xlab("Annual investment ($ Million)")+ylab("Proportion of 1750 habitat")+xlim(c(0,20))



curves.Vvals=data.frame(apply(curves.Vvals,2,function(x){ifelse(is.finite(x),x,1)}))

perfplot.Vvals=data.frame(curves$Cost,melt(curves.Vvals))
names(perfplot.Vvals)=c("Cost","Species","Value")

pp.Vvals=ggplot(data=perfplot.Vvals)+geom_line(aes(y=Value,x=Cost,group=Species,col=Species))
pp.Vvals+ggtitle("Habitat Value")+ theme(legend.position="none")

plants=rownames(features)[features$group=="plants"]
birds=rownames(features)[features$group=="birds"]
mammals=rownames(features)[features$group=="mammals"]


pp.Vvals.plants=ggplot(data=perfplot.Vvals[perfplot.Vvals$Species%in%plants,])+geom_line(aes(y=Value,x=Cost,group=Species,col=Species))
pp.Vvals.plants+ggtitle("Habitat Value: VROT plants")+ theme(legend.position="none")

pp.Vvals.birds=ggplot(data=perfplot.Vvals[perfplot.Vvals$Species%in%birds,])+geom_line(aes(y=Value,x=Cost,group=Species,col=Species))
pp.Vvals.birds+ggtitle("Habitat Value: All birds")+ theme(legend.position="none")

pp.Vvals.mammals=ggplot(data=perfplot.Vvals[perfplot.Vvals$Species%in%mammals,])+geom_line(aes(y=Value,x=Cost,group=Species,col=Species))
pp.Vvals.mammals+ggtitle("Habitat Value: All mammals")+ theme(legend.position="none")

totarea=sum(values(clip),na.rm=T)*225*225/10000

extinct.probs=1-curves.Vvals
perfplot.ext=data.frame(curves$Cost,(1-curves[,1])*totarea,apply(extinct.probs,1,sum))
names(perfplot.ext)=c("Cost","Area","Extinctions")
perfplot.ext$Cost=perfplot.ext$Cost/1000000
pp.ext=ggplot(data=perfplot.ext)+geom_line(aes(y=Extinctions,x=Cost),lwd=2)
pp.ext+ggtitle("Expected number of extinctions")+ theme_bw(base_size=16)+xlab("Annual investment ($ Million)")+ylab("Expected number of extinctions (50 years)")

perfplot.ext$Extinctions.avoided=max(perfplot.ext$Extinctions)-perfplot.ext$Extinctions

pp.ext=ggplot(data=perfplot.ext)+geom_line(aes(y=Extinctions.avoided,x=Cost),lwd=2)
pp.ext+ggtitle("Number of avoided extinctions")+ theme_bw(base_size=16)+xlab("Annual investment ($ Million)")+ylab("Avoided extinctions (50 years)")


extant.probs=curves.Vvals
perfplot.ext=data.frame(curves$Cost,(1-curves[,1]),apply(extant.probs,1,sum)/ncol(extant.probs))
names(perfplot.ext)=c("Cost","Area","Species")
perfplot.ext$Cost=perfplot.ext$Cost/1000000
pp.ext=ggplot(data=perfplot.ext)+geom_line(aes(y=Species,x=Cost),lwd=2)
pp.ext+ggtitle("Expected proportion of species persisting for 50 years")+ theme_bw(base_size=16)+xlab("Annual investment ($ Million)")+ylab("Proportion of extant species")

pp3=pp.ext+ theme_bw(base_size=16)+xlab("Annual investment ($ Million)")+ylab("Relative Biodiversity Value")

pp.extVarea=ggplot(data=perfplot.ext)+geom_line(aes(y=Species,x=Area),lwd=2)
pp6=pp.extVarea+ggtitle("Biodiversity Value vs Area treated")+ theme_bw(base_size=16)+xlab("Total area treated (sq. km)")+ylab("Relative Biodiversity Value")
pp6=pp.extVarea+ theme_bw(base_size=16)+xlab("Proportion of spatial priorities")+ylab("Relative Biodiversity Value")


costs.z=data.frame(values(zrank),values(bestact),matrix(0,ncell(zrank),ncol=nlayers(costs)),costvals)
names(costs.z)[1:2]=c("zrank","bestact")
costs.z=costs.z[-which(is.na(costs.z$zrank)),]

for(ln in colnames(totben.vals)){
  colnum=match(ln,colnames(totben.vals))
  costs.z[,ln]=ifelse(costs.z$bestact==colnum+1,costs.z$costvals,0)
}

costs.z=costs.z[order(costs.z$zrank,decreasing=T),]

dumcounts=table(costs.z$bestact)
dumcounts=dumcounts-dumcounts

seq=round(seq(0.05,1,by=0.05)*nrow(costs.z))

cumcounts=dumcounts
cumcosts=rep(0,ncol(costs.z)-2)
for(i in seq){
  newrow=dumcounts
  newvals=table(costs.z$bestact[1:i])
  newrow[names(newvals)]=newvals
  cumcounts=rbind(cumcounts,newrow)
  cv=costs.z[1:i,3:ncol(costs.z)]
  cumcosts=rbind(cumcosts,apply(cv,2,sum,na.rm=T))
}



cumprops=cumcounts/apply(cumcounts,1,sum)

cumcounts=cumcounts[,-1]
cumprops=cumprops[,-1]
colnames(cumprops)=colnames(cumcounts)=names(costs.all)[as.numeric(colnames(cumprops))]
head(costs.z)

cumcounts=cumcounts*225*225/1000000  #convert to square kms
cumcosts=cumcosts/1000000 #convert to millions of dollars
cumpdat=data.frame(c(0,seq),cumcosts[,"costvals"],melt(cumcosts[,colnames(cumcounts)])[,-1],melt(cumcounts)[,3])
#cumpdat=data.frame(c(0,seq),cumcosts[,"costvals"],cumcosts[,colnames(cumcounts)],melt(cumcounts))
names(cumpdat)=c("seq","cost","action","cost.a","area")
cumpdat$seq=cumpdat$seq/nrow(costs.z)

omitlevs=1#c(1,2)

pp1=ggplot(cumpdat,aes(x=cost,y=area,fill=action))+geom_area()+xlab("Annual Investment ($ Million)")+ylab("Area Treated (x1000 sq.km)")+ theme_bw(base_size=16)+theme(legend.title=element_blank(),legend.position=c(.2,.8))+
  scale_fill_manual(values = mycols[levels(bestact)[[1]]$ID][-omitlevs],labels=levs$action[-omitlevs]) 

pp1.m=pp1+theme(legend.position="none")


pp4=ggplot(cumpdat,aes(x=seq,y=area,fill=action))+geom_area()+xlab("Proportion of spatial priorities")+ylab("Area Treated (sq.km)")+ theme_bw(base_size=16)+theme(legend.title=element_blank(),legend.position=c(.2,.8))+
  scale_fill_manual(values = mycols[levels(bestact)[[1]]$ID][-omitlevs],labels=levs$action[-omitlevs]) 

pp4.m=pp4+theme(legend.position="none")



pp5=ggplot(cumpdat,aes(x=seq,y=cost.a,fill=action))+geom_area()+xlab("Proportion of spatial priorities")+ylab("Annual Investment ($ Million)")+ theme_bw(base_size=16)+theme(legend.title=element_blank(),legend.position=c(.3,.55))+
  scale_fill_manual(values = mycols[levels(bestact)[[1]]$ID][-omitlevs],labels=levs$action[-omitlevs]) 

pp5.m=pp5#+theme(legend.position="none")


pp2=ggplot(cumpdat,aes(x=cost,y=cost.a,fill=action))+geom_area()+xlab("Annual Investment ($ Million)")+ylab("Annual Investment ($ Million)")+ theme_bw(base_size=16)+theme(legend.title=element_blank(),legend.position=c(.2,.8))+
  scale_fill_manual(values = mycols[levels(bestact)[[1]]$ID][-omitlevs],labels=levs$action[-omitlevs]) 

pp2.m=pp2+theme(legend.position="none")


multiplot(pp1.m,pp2.m,pp3,pp4.m,pp5.m+theme(legend.position="none"),pp6,cols=2)

multiplot(pp1.m,pp2.m,pp3,pp4.m,pp5.m+theme(legend.position="none"),pp6,cols=2)
multiplot(pp1.m,pp2.m,pp3,cols=1)

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend <- g_legend(pp1+theme(legend.position=c(0.5,0.5))) 
grid.draw(legend) 

save.image(paste0("S:/SAN_Projects/SMP/Demo_Sept2015/SMP_Demo_Results_",Sys.Date(),".RData"))

### some totben and BCratio plots...
par(mar=c(2,2,2,2))
plot(totbens.tr[[8]],axes=F,legend=F,main="Potential Benefits: Fox and Cat Control",box=F)

plot(totbens.tr[[4]],axes=F,legend=F,main="Potential Benefits: Rabbit Warren Control",box=F)
plot(BCratio.tr[[4]],axes=F,legend=F,main="Benefit/Cost Ratio: Rabbit Warren Control",box=F)

plot(totbens.tr[[9]],axes=F,legend=F,main="Potential Benefits: Manage grazing and weeds",box=F)
plot(BCratio.tr[[9]],axes=F,legend=F,main="Benefit/Cost Ratio: Manage grazing and weeds",box=F)

plot(totbens.tr[[11]],axes=F,legend=F,main="Potential Benefits: Revegetation",box=F)
plot(BCratio.tr[[11]],axes=F,legend=F,main="Benefit/Cost Ratio: Revegetation",box=F)


## species specific outputs

egnames=c("Bush_Stone_curlew_Spp10174","Brush_tailed_Phascogale_Spp11017","Long_nosed_Bandicoot_Spp11097")
egnames=c("Turquoise_Parrot_Spp10302","Smoky_Mouse_Spp11458","Cape_Barren_Goose_Spp10198")

egnames=c("Bush_Stone_curlew_Spp10174","Cape_Barren_Goose_Spp10198")

for(spname in egnames[1:3]){
  setwd(rasterStack.folder)
  stacknames=list.files(pattern="_deltas")
  
  load(paste0(spname,"_deltas.Rstack"))
  #plot(deltas)
  names(deltas)=gsub("reveg","Reveg",names(deltas))
  bens=values(deltas)
  bcratio=values(deltas)/values(costs)[,names(deltas)]
  
  for(runit in 1:2){
    usevals=bcratio[viccells,]
    if(runit==2){usevals=bens[viccells,]}
    usevals[,"Clearing"]=usevals[,"Harvesting"]=usevals[,"Fuelred"]=0
    argmax=function(x){
      maxx=max(x,na.rm=T)
      ifelse(maxx>0,which(x==maxx)[1],0)
    }  
    maxs.spp=apply(usevals,1,argmax)
    #maxs.spp[apply(ifelse(is.na(usevals),0,usevals),1,sum)==0]=0
    bestact.spp=clip
    values(bestact.spp)[viccells]=maxs.spp+1
    bestact.spp=ratify(bestact.spp)
    levs.spp=levels(bestact.spp)[[1]]
    levs.spp$action=colnames(usevals)[levs.spp$ID]
    levs.spp$action=paste("Manage",levs.spp$action)
    levs.spp$action=gsub("Manage Reveg","Revegetate",levs.spp$action)
    levs.spp$action=gsub("Fox.Cats","Foxes & Cats",levs.spp$action)
    levs.spp$action=gsub("Manage Clearing","Control Clearing",levs.spp$action)
    levs.spp$action=gsub("Manage Dom.Grazing","Remove Stock",levs.spp$action)
    levs.spp$action=gsub("Manage Graz.Weeds","Manage Stock and  Weeds",levs.spp$action)
    levs.spp$action=gsub("Manage None","Urban area",levs.spp$action)
    levs.spp$action=gsub("Manage All","Manage all threats",levs.spp$action)
    levs.spp$action=gsub("Manage Fuelred","Burn below min TFI",levs.spp$action)
    levels(bestact.spp)=levs.spp
    
    delta2=addLayer(None,deltas)
    names(delta2)[1]="None"
    costs.best.spp=stackSelect.fast(costs.all[[names(delta2)]],values(bestact.spp))
    
    maxben.spp=stackSelect.fast(delta2,values(bestact.spp))
    #values(bestact.spp)=ifelse(values(bestact.spp)==1,NA,values(bestact.spp))
    #values(bestact.spp)[1:nrow(levels(bestact)[[1]])]=levels(bestact)[[1]]$ID
    pname1=c("Most cost-efficient actions","Most beneficial actions")[runit]
    plot(bestact.spp,breaks=c(levels(bestact.spp)[[1]]$ID-0.5,20),col=mycols[levels(bestact.spp)[[1]]$ID],legend=F,axes=F,main=pname1)
    legend(x="topright",legend=levels(bestact.spp)[[1]]$action,fill=mycols[levels(bestact.spp)[[1]]$ID])
    writeRaster(bestact.spp,file=file.path("S:/SAN_Projects/SMP/Demo_Oct2015/geotiffs",paste0(spname,"_",pname1,".tif")),format="GTiff",overwrite=T)
    
    spq=quantile(values(maxben.spp),0.975,na.rm=T)
    maxben.spp.trunc=maxben.spp
    values(maxben.spp.trunc)=ifelse(values(maxben.spp)>spq,spq,values(maxben.spp))
    pname2=c("Benefit of cost-efficient actions","Maximum benefit")[runit]
    plot(maxben.spp.trunc,axes=F,legend=F,main=pname2)
    writeRaster(maxben.spp.trunc,file=file.path("S:/SAN_Projects/SMP/Demo_Oct2015/geotiffs",paste0(spname,"_",pname2,".tif")),format="GTiff",overwrite=T)
    
    
    if(runit==1){
      bcrat=maxben.spp/costs.best.spp
      rank.spp=rank(-values(bcrat))
      is.fin=which(is.finite(values(maxben.spp)))
      maxben.cum=cumsum(values(maxben.spp)[is.fin][order(-values(bcrat)[is.fin])])
      maxben.cum=(maxben.cum+noact.sums[spname])/base.sums[spname]
      cost.cum=cumsum(values(costs.best.spp)[is.fin][order(-values(bcrat)[is.fin])])
      cost.cum=cost.cum/1000000
    }
    if(runit==2){
      best.cum=cumsum(sort(values(maxben.spp),decreasing=T))
      best.cum=(best.cum+noact.sums[spname])/base.sums[spname]
      fincosts=which(is.finite(values(costs.best.spp)))
      best.cost.cum=cumsum(values(costs.best.spp)[fincosts][order(-values(maxben.spp)[fincosts])])/1000000
      
      seq1000=seq(1,length(maxben.cum),by=1000)
      win.metafile(file.path("S:/SAN_Projects/SMP/Demo_Oct2015/geotiffs",paste0(spname,"_cumulativeBenefits.wmf")))
      costvals=curves$Cost/1000000
      maxcost=max(max(costvals,na.rm=T),min(max(cost.cum,na.rm=T),max(best.cost.cum,na.rm=T)))
      plot(best.cum[seq1000]~best.cost.cum[seq1000],type="l",lwd=2,main=paste("Cumulative Benefit"),xlab="Annual Investment ($ Million)",ylab="Proportion of 1750 distribution",ylim=c(0,1))
      lines(curves.Rvals[,spname]~costvals,col="blue",lwd=2)
      lines(maxben.cum[seq1000]~cost.cum[seq1000],col="red",lwd=2)
      dev.off()
    }
  }
}
### species over thresholds...

thresh=0.5
nt25=apply(ifelse(curves.Rvals>.25,1,0),1,sum)#/ncol(curves.Rvals)
nt50=apply(ifelse(curves.Rvals>.5,1,0),1,sum)#/ncol(curves.Rvals)
nt75=apply(ifelse(curves.Rvals>.75,1,0),1,sum)#/ncol(curves.Rvals)
costvals=curves$Cost/1000000

nt50a=nt50-min(nt50)
nt25a=nt25-min(nt25)
nt75a=nt75-min(nt75)

plot(nt50a~costvals,type="l",lwd=2,xlab="Annual Investment ($ Million)",ylab=paste("Number of additional species exceeding threshold"))
lines(nt25a~costvals,lwd=2,col="red")
lines(nt75a~costvals,lwd=2,col="blue")
text(500,150,"25% of 1750 distribution",col="red")
text(500,180,"75% of 1750 distribution",col="blue")
text(500,205,"50% of 1750 distribution")

plot(nt50a~costvals,xlim=c(0,200),type="l",lwd=2,xlab="Annual Investment ($ Million)",ylab=paste("Number of additional species exceeding threshold"))
lines(nt25a~costvals,lwd=2,col="red")
lines(nt75a~costvals,lwd=2,col="blue")
text(150,130,"25% of 1750 distribution",col="red")
text(150,165,"75% of 1750 distribution",col="blue")
text(150,205,"50% of 1750 distribution")


thresh=0.1
nt=apply(ifelse(curves.deltas.prop>thresh,1,0),1,sum)/ncol(curves.Rvals)
plot(nt~curves$Cost)#,ylim=c(min(nt),1))


cats=c(0,10,25,50,100)

costmill=curves$Cost/1000000
cinds=NULL
for(c in cats){
  diffs=abs(costmill-c)
  cinds=c(cinds,which(diffs==min(diffs,na.rm=T)))}

cdat=curves.Rvals[cinds,]
#cdat=curves.Rvals.relnow[cinds,]
library(reshape2)
cdat.df=data.frame(melt(cdat))
cdat.df$Investment=rep(cats,ncol(cdat))

nowrows=data.frame(names(current.sums),current.sums/base.sums[names(current.sums)],"now")
names(nowrows)=names(cdat.df)

#cdat.df=rbind(cdat.df,nowrows)

bp=ggplot(cdat.df,aes(factor(Investment),value))+geom_violin(fill="grey")
#bp=ggplot(cdat.df,aes(factor(Investment),value))+geom_violin(fill="grey")
bp2=bp+ylab("Net Outcome (prop. of 1750 distribution)")+xlab("Annual Investment ($ Million)")



cdat.delta=curves.deltas.prop[cinds,]
#cdat.delta=curves.deltas.prop.relnow[cinds,]
library(reshape2)
cdat.delta.df=data.frame(melt(cdat.delta))
cdat.delta.df$Investment=rep(cats,ncol(cdat.delta))

bp.delta=ggplot(cdat.delta.df,aes(factor(Investment),value))+geom_boxplot(fill="grey")
bp.delta=ggplot(cdat.delta.df,aes(factor(Investment),value))+geom_violin(fill="grey")
bp1=bp.delta+ylab("Benefit (as prop. of 1750 distribution)")+xlab("Annual Investment ($ Million)")

#bp.delta=ggplot(cdat.delta.df,aes(factor(Investment),value))+geom_violin(fill="grey")

multiplot(bp1,bp2)


###make a total risk map (weighted sum of expect change in hab value with no action)
cl <- makeCluster(nnodes)
registerDoParallel(cl)

pb=tkProgressBar("Risk map",label="species",min=0,max=length(species),initial=0)
#riskmap=foreach(spp=iter(species)[1:3],.packages="raster",.combine=sum)%do%{
riskmap=rep(0,ncell(clip))
for(spp in species){
  stackname=file.path("C:/SMPDemo/Demo/RasterStacks/",paste0(spp,"_habvals.Rstack"))
  load(stackname)
  risk=habvals[["None"]]-baseHDM[[spp]]
  
  rweight=(sum(values(baseHDM[[spp]]),na.rm=T)/base.sums[spp])^-.25
  rweight=ifelse(is.finite(rweight),rweight,1)
  newriskvals=na0(values(risk))*rweight
  if(!is.na(sum(newriskvals))){
    riskmap=riskmap+newriskvals}
  if(is.na(sum(riskmap))){break()}
  setTkProgressBar(pb,match(spp,species))
}

stopCluster(cl)

risk.map=clip
values(risk.map)=riskmap
risk.map=risk.map*clip
plot(risk.map)
writeRaster(risk.map,file=file.path(output.folder,paste0("weighted.risk.map_",date,".tif")),format="GTiff")



####plot truncated totben and BC ration maps as png's

for(lname in names(BCratio.tr)){
  png(file.path("S:/SAN_Projects/SMP/Demo_Oct2015/Benefit_BCratio_pngs",paste0(lname,"_BCratio.png")),width=1100,height=650)
  par(mar=c(1,1,1,1))
  plot(BCratio.tr[[lname]],axes=F,legend=F,box=F)
  dev.off()
  png(file.path("S:/SAN_Projects/SMP/Demo_Oct2015/Benefit_BCratio_pngs",paste0(lname,"_benefits.png")),width=1100,height=650)
  par(mar=c(1,1,1,1))
  plot(totbens.tr[[lname]],axes=F,legend=F,box=F)
  dev.off()
  if(lname%in%names(Threat.maps)){
    png(file.path("S:/SAN_Projects/SMP/Demo_Oct2015/Benefit_BCratio_pngs",paste0(lname,"_Threat.png")),width=1100,height=650)
    par(mar=c(1,1,1,1))
    plot(Threat.maps[[lname]],axes=F,legend=F,box=F)
    dev.off()
  }
  png(file.path("S:/SAN_Projects/SMP/Demo_Oct2015/Benefit_BCratio_pngs",paste0(lname,"_Costs.png")),width=1100,height=650)
  par(mar=c(1,1,1,1))
  plot(costs[[lname]],axes=F,legend=F,box=F)
  dev.off()
}

current.sums=colSums(values(baseHDM)[viccells,],na.rm=T)

}

