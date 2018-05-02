

first.start=Sys.time()

ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw(base_size=16) }


library(raster)
library(rasterVis)

library(ggplot2)

library(openxlsx)
library(doParallel)
library(tcltk)

source('O:/SAN_Projects/SMP/AAA_Rscripts/SMP_Rfunctions.R')
source('O:/SAN_Projects/SMP/AAA_Rscripts/SMP_RUN_functions.R')
rasterOptions(tmpdir="C:/temp/")
#rasterOptions(tmpdir="E:TEMP/")
input.file="O:/SAN_Work/SMPv1.0/runfiles/SMP_input files_sei.xlsx"

memory.limit(100000000)


## SET INPUTS and analysis settings


settings=read.xlsx(input.file,sheet="settings",rowNames=T)
nnodes=20#as.numeric(settings["MCnodes","Value"])

threats.info=read.xlsx(input.file,sheet="threats",rowNames=T)
threats=rownames(threats.info)

action.info=read.xlsx(input.file,sheet="actions",rowNames=T)
all.actions=rownames(action.info)
actions=all.actions#[-grep("Reveg",all.actions)]

species.info=read.xlsx(input.file,sheet="species",rowNames=T)

#load vic and NV mask rasters..

vic=raster("E:/SMPv1.2/Inputs/vicmask225vg.tif")
vicnv=raster("E:/SMPv1.2/Inputs/nv225vg.tif") # or use vicNV225_prop.tif

clip=vic
clipvals=values(clip)
viccells=which(clipvals==1)
nvcells=which(values(vicnv)==1) 

output.folder=settings["output.folder","Value"]
dir.create(output.folder)
Zwd=settings["Zonation.directory","Value"]
dir.create(Zwd)
dir.create(paste0(Zwd,"/inputs"))
#load species HDMs
#setwd(settings["SDM.folder","Value"])
#baseHDM=stack(rownames(species.info))
#baseHDM=raster.multiply(baseHDM,0.01*values(clip))
#save(baseHDM,file=file.path(output.folder,"baseHDM.Rstack"))
#load("baseHDM.Rstack")

species.full=rownames(species.info)
species=species.info$taxoncode

ncells=ncell(vic)


##load threat maps (as stack)
Threat.maps=stack(stack(threats.info[,1]))
names(Threat.maps)=threats #rownames(threats.info)

rasterStack.folder=paste0(output.folder,"/RasterStacks")
dir.create(output.folder)
dir.create(rasterStack.folder)
setwd(rasterStack.folder)

AddSpecies=settings["AddSpecies","Value"]

do.species=species

#load cost maps (as stack)
costs=stack(action.info$CostFile)
names(costs)=all.actions# rownames(action.info)
for(cname in names(costs)){
  if(cname%in%names(Threat.maps)){
    values(costs[[cname]])=ifelse(values(costs[[cname]])<0.02,0.02,values(costs[[cname]]))
    values(costs[[cname]])=ifelse(values(Threat.maps[[cname]])>0,values(costs[[cname]]),NA)
  }
}
grazweedcost=ifelse(values(Threat.maps[["Weeds"]])>0.1,ifelse(values(Threat.maps[["Dom.Grazing"]])>0.5,values(costs[["Graz.Weeds"]]),values(costs[["Weeds"]])),values(Threat.maps[["Dom.Grazing"]])*values(costs[["Dom.Grazing"]]))
values(costs[["Graz.Weeds"]])=ifelse(values(Threat.maps[["Weeds"]])>0.1 & values(Threat.maps[["Dom.Grazing"]])>0.5,values(costs[["Graz.Weeds"]]),NA)
#values(costs[["All"]])=na0(values(costs[["Fox.Cats"]]))+na0(grazweedcost)+rowSums(values(costs[[names(costs)[names(costs)%in%c("Foxes","Cats","Dom.Grazing","Weeds","Reveg","All","Harvesting")==FALSE]]]),na.rm=T)
#values(costs[["All"]])=values(costs[["All"]])*values(clip)
mincells=action.info$mincells
names(mincells)=rownames(action.info)

##addd feral costs to reveg
#costs[["Reveg"]]=costs[["Reveg"]]+0.5*(costs[["Weeds"]]+0.5*costs[["Fox.Cats"]])

#set zero cost values to NA...
ztoNA=function(x){ifelse(x<=0,NA,x)}
values(costs)=apply(values(costs),2,ztoNA)

for(actn in names(which(mincells>1))){
  nbv=round(sqrt(mincells[actn]))
  if(nbv%%2==0){nbv=nbv+1}
  costs[[actn]]=focal(costs[[actn]],w=matrix(1,nbv,nbv),fun=mean,na.rm=T,pad=T)*clip
}

## dividing predator control values by 10....(temporary fix!!)
# costs[["Foxes"]][Threat.maps[["Foxes"]]==0]=NA
# costs[["Fox.Cats"]][Threat.maps[["Foxes"]]==0]=NA
# 
# costs[["Foxes"]]=costs[["Foxes"]]/10
# costs[["Cats"]]=costs[["Cats"]]/10
# costs[["Fox.Cats"]]=costs[["Fox.Cats"]]/10

## making clearing costs = average TfN covenent cost, and excluding already protected land
#costs[["Clearing"]][costs[["Clearing"]]]=5.0625*750/50 # TfN covenent per ha spread over 50 years (X pixel size)
costs[["Clearing"]][Threat.maps[["Clearing"]]==0]=NA

#make costs per 5.0625 ha pixel
for(lyr in names(costs)[-1]){costs[[lyr]]=costs[[lyr]]*5.0625/4} #Matts surfaces were cost per 4 ha pixel
costs[["Clearing"]]=costs[["Clearing"]]*5.0625

writeRaster(costs,file=file.path(output.folder,"SMPCosts.tif"),format="GTiff",overwrite=T)
## calculate local benefits of each action set..
# Scale by max. possible benefit - based on benefit fuction..

# calc max possibile benefit and smoothed delta layers

features=species.info
features$sum.delta.max=features$maxben=features$minBF=features$minprop=features$maxprop=features$sumreveg=NA
#convert 1750 sum from ha to # of pixels.
features$sum1750cells=as.numeric(features$sum1750ha)/(prod(res(clip))/10000)

dudhabvals=NULL
reveg0=clip-1

## Create 'delta' staks for each species, and calculated weighted sum of delta values (benefits)

## This is the basis of stage 1 prioritization. Nspecies X Nactions rasters built- long processing times!!
BenefitDir=settings["BenefitRastersDirectory","Value"]

colbens=collateBenefits(
  do.species=do.species,
  sppInfo=features,# list of taxoncodes of the species to include
  BenefitDir=BenefitDir, # directory where benefit rasters are stored
  rasterStack.folder=rasterStack.folder, #folder where (smoothed) benefit raster stacks will be stored (will be created if doesn't exist)
  actions=actions,#action names
  nvactions=actions[-grep("Reveg",actions)],# actions that apply to native veg (generally all but reveg)
  AddSpecies=F,addTo=NULL,#"E:/SMPv1.0/Results/sensitivity/RasterStacks/totben.vals_3.Rmat",  #TRUE if adding species to existing totben.vals table (addTo = full path to that R object)
  nnodes=nnodes, #number of parallel nodes to use...
  mask=raster("O:/SAN_Projects/SMP/vicmask225.tif"), # analysis mask (cells with value >0 included in analysis)
  vegmask=raster("O:/SAN_Projects/SMP/vicNV225_binary.tif"), # native veg mask (> 0 for native veg)
  costs=costs,
  smoothbenefits=F, # smooth benefits over minimum area?
  mincells=mincells,# names vector of min. no. of pixels over which an action should be performed.
  onParError="remove",maxgpsize=400
  
)

##append 'features' df with weights etc for later use
#load the weightings file
load(colbens$zweights.file)
gotem=rownames(weightings)[rownames(weightings)%in%features$taxoncode]
features[match(as.numeric(gotem),features$taxoncode),names(weightings)]=weightings[gotem,]
#save(smoothit[[2]],file=paste0("smoothit_",gp,".R"))

minprop.fill=ifelse(is.na(features$minprop),mean(features$minprop,na.rm=T),features$minprop)
maxprop.fill=ifelse(is.na(features$maxprop),mean(features$maxprop,na.rm=T),features$maxprop)

bfparams=bf.adj.multi(features$Vx,minprop.fill,maxprop.fill)
features=data.frame(features,bfparams)
features$Zweight=features$weight*features$weight.1
save(features,file=file.path(output.folder,"features.Rmat"))


###load saved totben.vals and features.Rmat
load(colbens$totben.file)
load(file.path(output.folder,"features.Rmat"))

## create raster stacks of summed benefits and BC ratios..
totbens=clip-1
for(a in 2:length(actions)){totbens=addLayer(totbens,(clip-1))}
names(totbens)=actions
values(totbens)[viccells,]=totben.vals
BCratio.all=totbens
values(BCratio.all)[viccells,]=totben.vals/values(costs)[viccells,colnames(totben.vals)]

writeRaster(totbens,file=file.path(output.folder,"totalbenefits.tif"),format="GTiff",overwrite=T)
writeRaster(BCratio.all,file=file.path(output.folder,"BCratios.tif"),format="GTiff",overwrite=T)
tben.names=names(totbens)
save(tben.names,file=file.path(output.folder,"BCratios.tifnames.R"))

save(totbens,file=file.path(output.folder,"totbens.Rstack"))
save(BCratio.all,file=file.path(output.folder,"BCratio.Rstack"))

actions.scen=actions[actions%in%c("Harvesting","Phytoph","Fuelred")==F]
scenario.name=settings["scenario.name","Value"]
#subset actions for scenario...

totben.vals.scen=totben.vals[,actions.scen]

#make benefits of dom.grazing control zero on public land (can't stop grazing on public land...!?)
pub225=raster("S:/SAN_Work/SMPv1.0/Inputs/publicLand225.tif")
pub225=resample(pub225,clip)
pubvic=values(pub225)[viccells]
#BCvals[,"Dom.Grazing"]=0

totben.vals.scen[pubvic==1,"Dom.Grazing"]=0
totben.vals.scen[pubvic==1,"Graz.Weeds"]=0

costs.scen=costs[[actions.scen]]

cost.vals=values(costs.scen)[viccells,actions.scen]

BCvals=totben.vals.scen/cost.vals


#find best local action
#final step of state 1 prioritization

clearing.probs=values(Threat.maps[["Clearing"]])[viccells]
clearing.probs=1-((1-clearing.probs)^0.05)^50   #convert 20 year clearing prob to 50 year clearing prob


#get best actions...

nogo.combos=rbind(c("Weeds","Dom.Grazing","Graz.Weeds"),c("Cats","Foxes","Fox.Cats"))  #makes some combinations of actions out of bounds (make user friendly in app?)


bestactions.scen=getBestActions(benefits=totben.vals.scen,
                                costs=cost.vals,
                                clearing.probs=clearing.probs,
                                nogo.combos=nogo.combos,
                                clearing.ind=1,all=999)

load(bestactions.scen$raster)
lookup=bestactions.scen$lookup

## make top combinations map....




date=Sys.Date()
writeRaster(bestact,file=paste0(output.folder,"/SMPbestaction",date,".tif"),format="GTiff",overwrite=T)


##END STAGE 1 PRIORITIZATION

#copy best action cost tif to Zwd folder

file.copy(bestactions.scen$cost.raster,file.path(Zwd,"/inputs/SMP_costbest.tif"),overwrite=T)


#make best action benefits rasters...


weight.rows=getBAbenefits(bestact,lookup,rasterStack.folder,features,species=do.species,BenefitDir,actions=actions.scen,
                          species.done=NULL,AddSpecies=FALSE,Zwd=Zwd,mask=clip,nnodes=nnodes,outpath=output.folder,onParError="remove")







#Write Zonation run & setting files (incl. cell removel rule, warp factor, edge removal etc).

Zfiles=makeZfiles(Zwd,features,weight.rows,
    removal.rule=settings["removal.rule","Value"],
    warp=settings["warp","Value"],
    edge.removal=settings["edge.removal","Value"],
    edge.points=settings["edge.points","Value"]
)
#write Zonation .dat file

### RUN ZONATION.....

system(Zfiles$batfile)

#Re Run with reveg limit
Zfiles.reveg=Zreveglim(Zwd,Zfiles,reveg.lim=2)

system(Zfiles.reveg$batfile)

######  
  
  
