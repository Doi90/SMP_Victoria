
library(mboost)
data("bodyfat", package = "TH.data")
mod <- mboost(DEXfat ~ btree(age) + bols(waistcirc) + bbs(hipcirc),data = bodyfat)

library(raster)
library(dplyr)
library(tcltk)
library(doParallel)
library(rgdal)
library(randomForest)
library(openxlsx)
library(verification)

library(reshape2)

na0=function(x){ifelse(is.na(as.numeric(x)),0,as.numeric(x))}
sumna0=function(x){sum(na0(x))}

#load expert answers:
responses=read.xlsx("S:/SAN_Projects/SMP/ExpertElicitation/Answers/EEresults_all_13October2015.xlsx",sheet=1,rowNames=T)
#load site data (threats)
sitedata=read.xlsx("S:/SAN_Projects/SMP/ExpertElicitation/Answers/EEresults_all_13October2015.xlsx",sheet=2,rowNames=T)
#load taxon info (mainly used to group species into taxonomic groups...)
taxainfo=read.xlsx("S:/SAN_Projects/SMP/ExpertElicitation/Answers/EEresults_all_13October2015.xlsx",sheet=3,rowNames=F)
taxainfo=unique(taxainfo)

#compile expert data for modelling 
##NB might be better to import...S:/SAN_Projects/SMP/ExpertElicitation/EEdataViewer/response.df.Rdf

expertnames=unique(names(responses)[grep("mean_",names(responses))])
expertnames=expertnames[-grep("refmean",expertnames)]

for(expert in expertnames){
  mean_benefit=as.numeric(responses[,expert])-as.numeric(responses[,paste0("ref",expert)])
  responses[paste0(expert,"_benefit")]=mean_benefit
}

rownames(taxainfo)=taxainfo$TAXON_ID

sitedata[,1:ncol(sitedata)]=apply(sitedata[,1:ncol(sitedata)],2,na0)

taxa=unique(responses$TaxonCode)

response.melt=melt(responses[-grep("comments",names(responses))],id.vars=names(responses)[1:11],na.rm=T)
response.melt=response.melt[is.finite(as.numeric(response.melt$value)),]

#
group=taxainfo$TAXON_TYPE[match(response.melt$TaxonCode,taxainfo$TAXON_ID)]
group[grep("birds",group)]="Birds"
group[grep("Waders",group)]="Birds"
group[group%in%c("Dicotyledons","Monocotyledons","Fern and allies","Conifers")]="Plants"

response.melt$group=group

#subset EE data...using only mammals from here on...

response.mammals=response.melt[response.melt$group=="Mammals",]

response.mammals$variable=factor(gsub("2","",response.mammals$variable))
response.mammals=unique(response.mammals)

response.mammals.all=response.mammals
#remove refmeans & lower uppper estimates (for now)..
response.mammals.refmean=response.mammals[grep("refmean_",as.vector(response.mammals$variable)),]
response.mammals.refmean=response.mammals.refmean[response.mammals.refmean$action=="Do Nothing",]
response.mammals=response.mammals[grep("_benefit",as.vector(response.mammals$variable)),]
response.mammals=response.mammals[response.mammals$action!="Do Nothing",]


response.mammals$variable=factor(gsub("mean_","",gsub("_benefit","",response.mammals$variable)))
response.mammals.refmean$variable=factor(gsub("refmean_","",response.mammals.refmean$variable))

##Get species trait data

mamm.traits=read.csv("S:/SAN_Projects/SMP/Trait Data/MammalTraits_for clustering.csv")

mamm.traits$max.hr=apply(mamm.traits[,11:16],1,max,na.rm=T)
mamm.traits$min.hr=apply(mamm.traits[,11:16],1,min,na.rm=T)
mamm.traits$min.hr=ifelse(is.finite(mamm.traits$min.hr),mamm.traits$min.hr,NA)
mamm.traits$max.hr=ifelse(is.finite(mamm.traits$max.hr),mamm.traits$max.hr,NA)

mamm.traits$msm=ifelse(is.na(mamm.traits$MaleSexMat.min),mamm.traits$MSM.max,mamm.traits$MaleSexMat.min)
mamm.traits$msm=ifelse(is.na(mamm.traits$FSM.min),mamm.traits$FSM.max,mamm.traits$FSM.min)
mamm.traits$litter.size=ifelse(is.finite(mamm.traits$Litter.size.average),mamm.traits$Litter.size.average,apply(mamm.traits[,c("litter.size.min","litter.size.max")],1,mean,na.rm=T))
mamm.traits$litter.size=ifelse(is.nan(mamm.traits$litter.size),NA,mamm.traits$litter.size)
mamm.traits$Mlong=apply(mamm.traits[,grep("MaleLongevity",names(mamm.traits))],1,max,na.rm=T)
mamm.traits$Flong=apply(mamm.traits[,grep("FemaleLongevity",names(mamm.traits))],1,max,na.rm=T)
mamm.traits$Mlong=ifelse(is.finite(mamm.traits$Mlong),mamm.traits$Mlong,NA)
mamm.traits$Flong=ifelse(is.finite(mamm.traits$Flong),mamm.traits$Flong,NA)

mamm.traits=mamm.traits[]

threats=sitedata[response.mammals$ref.scenario,]


trait.cols=c(4,6,7,9,10,30,32,37:40) # using these traits...


traits=mamm.traits[match(response.mammals$TaxonCode,mamm.traits$TaxonCode),trait.cols]

#get scale, strata, climate change settings from scenario title

scale=as.numeric(gsub(".*_","",gsub("ha_.*","",response.mammals$title)))
scen=gsub("_st.*","",response.mammals$ref.scenario)
strata=as.numeric(gsub("_sc.*","",gsub(".*_st","",response.mammals$ref.scenario)))
nocc=rep(0,length(strata))
nocc[grep("_noCC_",response.mammals$ref.scenario)]=1

expert.refscen=paste0(response.mammals$ref.scenario,response.mammals$variable)
expert.refscen.refmean=paste0(response.mammals.refmean$ref.scenario,response.mammals.refmean$variable)
refmean=as.numeric(response.mammals.refmean$value)[match(expert.refscen,expert.refscen.refmean)]

#make data frame for input into models, here y (response.mammals$value) is the expected benefit change from reference (NoAction) prediction to with action prediction
mdat=data.frame(y=as.numeric(response.mammals$value),action=response.mammals$action,refmean,tempscale=response.mammals$tempscale,scale,threats,strata,nocc,traits)

fmla=paste("y~btree(",paste(names(mdat)[2:ncol(mdat)],collapse=","),")")
fmla2=gsub("action","species, action",fmla)


fmla=as.formula(paste(fmla,"+ brandom(species) + brandom(expert)"))  #+ brandom(scen) 
fmla2=as.formula(paste(fmla2," + brandom(expert)")) # + brandom(scen)

mdat=data.frame(mdat,expert=response.mammals$variable,scen,species=factor(response.mammals$TaxonCode))
vv=as.vector(mdat$expert)
#mdat$expert=factor(ifelse(vv!="plent" & vv!="llums","other",vv))

###Save workspace
save.image("S:/SAN_Projects/SMP/ModelEEdata_mammals.RData")

#build models...

#First model (fmla) has random species effects...used to predict to species not included in oringal expert elicitation (based on traits only) 

mod=mboost(fmla,data=mdat,control=boost_control(mstop=1000,trace=T))

#get fitted values
preds=predict(mod,newdata=mdat,which=1)
preds=preds+attr(preds,"offset")
plot(preds~factor(mdat$y))

## use CV to find optimal boosting iteration...
myApply <- function(X, FUN, ...) {
  myFun <- function(...) {
    library("mboost") # load mboost on nodes
    FUN(...)
  }
  ## further set up steps as required
  parLapply(cl = cl, X, myFun, ...)
}
cl <- makeCluster(10)
optn=cvrisk(mod, papply = myApply)
stopCluster(cl)

mod[991]

save(mod,file="S:/SAN_Projects/SMP/MboostModels/mboost_mammals_benefits_randspecies.R")

#model 2 (fmla2) has species as fixed effect...used for predicting to species that have EE data

mod2=mboost(fmla2,data=mdat,control=boost_control(mstop=1000,trace=T))

cl <- makeCluster(10) # e.g. to run cvrisk on 25 nodes via PVM
cvrisk(mod2, papply = myApply)
stopCluster(cl)

mod2[421]

preds=predict(mod2,newdata=mdat,which=1)
preds=preds+attr(preds,"offset")
plot(preds~factor(mdat$y))

save(mod2,file="S:/SAN_Projects/SMP/MboostModels/mboost_mammals_benefits.R")

preds=predict(mod2,newdata=mdat,which=1,off2int=T)
preds=preds+attr(preds,"offset")
plot(preds~factor(mdat$y))



###NOW predict persistence probs conditional on No Action.

threats=sitedata[response.mammals.refmean$ref.scenario,]
traits=mamm.traits[match(response.mammals.refmean$TaxonCode,mamm.traits$TaxonCode),trait.cols]

scale=as.numeric(gsub(".*_","",gsub("ha_.*","",response.mammals.refmean$title)))
scen=gsub("_st.*","",response.mammals.refmean$ref.scenario)
strata=as.numeric(gsub("_sc.*","",gsub(".*_st","",response.mammals.refmean$ref.scenario)))
nocc=rep(0,length(strata))
nocc[grep("_noCC_",response.mammals.refmean$ref.scenario)]=1


mdat.na=data.frame(y=as.numeric(response.mammals.refmean$value),action=response.mammals.refmean$action,tempscale=response.mammals.refmean$tempscale,scale,threats,strata,nocc,traits)

fmla=paste("y~btree(",paste(names(mdat.na)[3:ncol(mdat.na)],collapse=","),")")
fmla2=gsub("tempscale","species, tempscale",fmla)
fmla=as.formula(paste(fmla,"+ brandom(species) + brandom(expert)"))  #+ brandom(scen) 
fmla2=as.formula(paste(fmla2," + brandom(expert)")) # + brandom(scen)

mdat.na=data.frame(mdat.na,expert=response.mammals.refmean$variable,scen,species=factor(response.mammals.refmean$TaxonCode))
vv=as.vector(mdat.na$expert)
#mdat.na$expert=factor(ifelse(vv!="plent" & vv!="llums","other",vv))

# random species model

mod.na=mboost(fmla,data=mdat.na,control=boost_control(mstop=1000,trace=T))

plot(fitted(mod.na)~factor(mdat.na$y))

preds=predict(mod.na,newdata=mdat.na,which=1)
preds=preds+attr(preds,"offset")
plot(preds~factor(mdat.na$y))

cl <- makeCluster(10) # e.g. to run cvrisk on 10 nodes via PVM
cvrisk(mod.na, papply = myApply)
stopCluster(cl)

mod.na[291]

save(mod.na,file="S:/SAN_Projects/SMP/MboostModels/mboost_mammals_NoAction_randspecies.R")

# fixed species model

mod2.na=mboost(fmla2,data=mdat.na,control=boost_control(mstop=1000,trace=T))

cl <- makeCluster(4) # e.g. to run cvrisk on 10 nodes via PVM
cvrisk(mod2.na, papply = myApply)
stopCluster(cl)

mod2.na[190]

plot(fitted(mod2.na)~factor(mdat.na$y))

save(mod2.na,file="S:/SAN_Projects/SMP/MboostModels/mboost_mammals_NoAction.R")

preds=predict(mod2.na,newdata=mdat.na,which=1,off2int=T)
preds=preds+attr(preds,"offset")
plot(preds~factor(mdat.na$y))


load("S:/SAN_Projects/SMP/MboostModels/mboost_mammals_benefits_randspecies.R")
load("S:/SAN_Projects/SMP/MboostModels/mboost_mammals_benefits.R")
load("S:/SAN_Projects/SMP/MboostModels/mboost_mammals_NoAction_randspecies.R")
load("S:/SAN_Projects/SMP/MboostModels/mboost_mammals_NoAction.R")
#map predictions

