
na0=function(x){ifelse(is.na(x),0,x)}
ztoNA=function(x){ifelse(x<=0,NA,x)}
##Zonation generalized benefit function (see Zonation manual)
BF=function(R,x=0.25,y=1,T=1,w1=1,w2=1){ifelse(R<=T,w1*((R/T)^x),w1+w2*(((R-T)/(1-T))^y))}

#additive benefit function...matrix based (no raster inputs)
AB=function(dvals,sp.weights,base.sums,noact.sums,costvals,expons){
  sumdelts=apply(dvals,2,sum,na.rm=T)
  sumdelts.noi=-sweep(dvals,2,sumdelts,"-")
  sumdelts.noi.noacts=sweep(sumdelts.noi,2,noact.sums,"+")
  exp.sums=(sumdelts+noact.sums)^expons
  exp.sums.noi=sweep(sumdelts.noi.noacts,2,expons,"^")
  diffs=-sweep(exp.sums.noi,2,exp.sums,"-")
  weight=sp.weights/(base.sums^expons)
  ivals=sweep(diffs,2,weight,"*")
  sum.ivals=apply(ivals,1,sum)
  sum.ivals.ce=sum.ivals/costvals
  remove=which(sum.ivals.ce==min(sum.ivals.ce,na.rm=T))[1]
  return(remove)
}

#additive benefit function, raster inputs, allows edge removal and warp factors

ABras=function(dvals.stack,sp.weights,base.sums,noact.sums,costvals,expons,edge.remove=FALSE,warp=1){
  cands=which(is.finite(values(dvals.stack[[1]])))
  if(edge.remove){
    edges=which(is.na(values(focal(dvals.stack[[1]],w=matrix(1,3,3),fun=sum))))
    cands=edges[edges%in%cands] 
  }
  dvals=values(dvals.stack)[cands,]
  sumdelts=colSums(dvals,na.rm=TRUE)
  sumdelts.noi=-sweep(dvals,2,sumdelts,"-")
  sumdelts.noi.noacts=sweep(sumdelts.noi,2,noact.sums,"+")
  exp.sums=(sumdelts+noact.sums)^expons
  exp.sums.noi=sweep(sumdelts.noi.noacts,2,expons,"^")
  diffs=-sweep(exp.sums.noi,2,exp.sums,"-")
  ivals=sweep(diffs,2,sp.weights,"*")
  sum.ivals=rowSums(ivals)
  sum.ivals.ce=sum.ivals/costvals[cands]
  remove=cands[order(sum.ivals.ce)[warp]]
  return(remove)
}

#additive benefit function, using Zonation's fast approximation. Raster based, with edge removal and warp factor

ABZ=function(dvals.stack,sp.weights,costvals,expons,edge.remove=FALSE,warp=1){
  cands=which(is.finite(values(dvals.stack[[1]])))
  if(edge.remove){
    edges=which(is.na(values(focal(dvals.stack[[1]],w=matrix(1,3,3),fun=sum))))
    cands=edges[edges%in%cands] 
  }
  dvals=values(dvals.stack)[cands,]
  Rvals=colSums(dvals,na.rm=TRUE)
  weight=sp.weights*expons*Rvals^(expons-1)
  ivals=sweep(dvals,2,weight,"*")
  sum.ivals=rowSums(ivals)
  sum.ivals.ce=sum.ivals/costvals[cands]
  remove=cands[order(sum.ivals.ce)[warp]]
  return(remove)
}

## set Zonation benefit function parameters so that mrginal loss of delta's approximates marginal loss for total habitat

bf.adj=function(exp,minprop,maxprop){
  params=c(0,0)
  pdiff=maxprop-minprop
  if(round(pdiff,6)!=0){
  seq=sort(seq(minprop,maxprop,length.out=10000))
  seq.adj=(seq-min(minprop,maxprop))/abs(pdiff)
  orig.slope=exp*(seq^(exp-1))
  params=lm(log(orig.slope[-1])~log(seq.adj[-1]))$coefficients
  }
  return(list(weight=abs(pdiff)*exp(params[1])/(params[2]+1),expon=(params[2]+1)))  
}


bf.adj.multi=function(exps,minprops,maxprops){
  res=NULL
  for(i in 1:length(exps)){
    res=rbind(res,data.frame(bf.adj(exps[i],minprops[i],maxprops[i])))
  }
  return(res)
}


### plot grids with labels...
plotgrid=function(grid,labels="vals",revy=TRUE,label.lu=NULL,truncscale=TRUE,alphavals=1,border="black",colregs=rev(terrain.colors(256))){
  myPanel <- function(x, y, z,...) {
    panel.levelplot(x,y,z,border=border,...)
    if(labels=="vals"){panel.text(x, y, ifelse(z==0,"",round(z,2)))}
    if(labels=="labels"){panel.text(x, y, label.lu[z],alpha=alphavals)}
  }
  x=c(1:ncol(grid))
  y=c(1:nrow(grid))
  gg=expand.grid(X=x,Y=y)
  if(revy){gg$Y=nrow(grid)-gg$Y+1}
  gg$Z=as.vector(grid)
  
  if(truncscale){colregs=colregs[1:ceiling(max(gg$Z)*256)]}
  levelplot(Z~X*Y,gg,panel=myPanel,colorkey=FALSE, col.regions=colregs,scales=list(draw=FALSE),xlab="",ylab="",alpha.regions=alphavals)
}

## Core area zonation...(check before using...not tested!)

CA=function(dvals,sp.weights,maxbens,costvals){
  dval.prop=sweep(dvals,2,apply(dvals,2,sum,na.rm=T),"/")
  weight=sp.weights*maxbens
  dval.prop.w=sweep(dval.prop,2,weight,"*")
  maxval=apply(dval.prop.w,1,max,na.rm=T)
  maxval.ce=maxval/costvals
  maxval.ce=ifelse(is.finite(maxval.ce),maxval.ce,NA)
  remove=which(maxval.ce==min(maxval.ce,na.rm=T))[1]
  return(remove)
}

stack.prod=function(X){
  ras=X[[1]]
  if(nlayers(X)>1){
  vals=values(ras)
  for(c in 2:nlayers(X)){vals=vals*values(X[[c]])}
  values(ras)=vals
  }
  return(ras)
}

subtract.raster=function(Stack,sublayer=names(Stack)[1]){
  newstack=Stack[[1]]
  values(newstack)=values(Stack[[2]])-values(Stack[[sublayer]])
  for(r in 3:nlayers(Stack)){
    newras=newstack[[1]]
    values(newras)=values(Stack[[r]])-values(Stack[[sublayer]])
    newstack=addLayer(newstack,newras)
  }
names(newstack)=names(Stack)[-1]
return(newstack)
}

raster.multiply=function(ras,c){for(l in 1:nlayers(ras)){values(ras[[l]])=values(ras[[l]])*c};return(ras)}


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

 
focal.mean=function(inraster){apply(inraster,3,mean,na.rm=T)}

majority=function(x){
  tab=table(na.omit(x))
  maj=names(which(tab==max(tab,na.rm=T)))
  if(is.numeric(x)){maj=as.numeric(maj)}
  return(maj)
}

max0=function(x){m=max(x,na.rm=T);m=ifelse(is.finite(m),m,0)}


#threat.impact=function(threat,vuln){exp(-vuln*threat)}
threat.impact=function(threat,vuln){threat*vuln + (1-threat)} # function giving proportional habitat value given "threat" level

map.impacts=function(threatmap,sp.vuln){
  vulras=threatmap
  values(vulras)=sp.vuln
  newprobs=overlay(threatmap,vulras, fun=threat.impact)
  return(newprobs)
}

sumRasVals=function(v1,v2){v1+v2}

stdz.stack=function(RasterStack,base.vals=NULL,npar.nodes=5,dir="E:/SMPDemo/SDMs/std/"){
  cl <- makeCluster(nnodes)
  registerDoParallel(cl)
  newStack=foreach(lyr=iter(names(RasterStack)),.packages="raster",.combine=sumRasVals)%dopar%
  {
  if(is.null(base.vals)){base.val=sum(values(RasterStack[[lyr]]),na.rm=T)}else{base.val=base.vals[lyr]}
  ras.std=RasterStack[[lyr]]/base.val
  save(ras.std,file=paste0(dir,lyr,".std.Rras"))
  vals=ifelse(is.na(values(ras.std)),0,values(ras.std))
  gc()
  vals
  }
  stopCluster(cl)
  return(newStack)
}

stdz.stack2=function(rasnames,base.std.sums,nzcount,npar.nodes=5,dir="E:/SMPDemo/SDMs/std/"){
  cl <- makeCluster(nnodes)
  registerDoParallel(cl)
  newStack=foreach(lyr=iter(rasnames),.packages="raster",.combine=c)%dopar%
{
  load(paste0(lyr,".std.Rras"))
  values(ras.std)=values(ras.std)/base.std.sums
  weight=sum(values(ras.std),na.rm=T)/nzcount[lyr]
  names(weight)=lyr
  rm(ras.std)
  gc()
  weight
}
stopCluster(cl)
return(newStack)
}

oneoff=function(){
newvals=base.std.sums
for(lyr in missingsums){
  load(paste0(lyr,".std.Rras"))
  newvals=newvals+ifelse(is.na(values(ras.std)),0,values(ras.std))
  rm(ras.std)
  gc()
}
}
#plot raster stack with common color gradient

plotStack=function(stack,ncols=255,col=rev(terrain.colors(ncols)),legend=FALSE){
  mybreaks=seq(min(minValue(stack)),max(maxValue(stack)),length.out=ncols)
  plot(stack,breaks=mybreaks,col=col,legend=legend)
}

###draw multi-panel plot from ggplot objects
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL,gridlayout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if(is.null(gridlayout)){gridlayout=grid.layout(nrow(layout), ncol(layout))}
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = gridlayout))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



addColorTable <- function(inRstName, outRstName, rat.df){
  library(rgdal)
  r<- readGDAL(inRstName)
  rat.df$color<- as.character(rat.df$color)
  rat.df$attribute<- as.character(rat.df$attribute)
  outRst <- writeGDAL(r, outRstName, type="Byte", 
                      colorTable=list(rat.df$color), 
                      catNames=list(rat.df$attribute), mvFlag=11L)
  return(raster(outRst))
}

#check how this adds cost to covenants
argmax.comb=function(ben,cost, clearing.prob=0,clearing.ind=1,all=999,nogos=nogo.combos){
  bestcol=bestid=tben=tcost=NA
  cands=which(ben/cost>0)
  pext=rep(1-clearing.prob,length(ben))
  pext[clearing.ind]=1
  if(length(cands)>0){
  ben=ben[cands]
  cost=cost[cands]
  pext=pext[cands]
  ranks=rank(cost/(ben*pext),ties.method="first")
  best=which(ranks==1)
  best.val=((ben*pext)/cost)[best]
  if(cands[best]==clearing.ind){pext=1}
  if(cands[best]!=all){
  if(all%in%cands){ben[cands==all]=0}
  for(cr in 1:nrow(nogos)){
  if(sum(cands[best]%in%nogos[cr,])>0){ben[cands%in%nogos[cr,]]=0}}
  ben.rem=ben
  cost.rem=0.5*cost
  cost.rem[which(cands==clearing.ind)]=cost[which(cands==clearing.ind)]
  ben.rem[best]=cost.rem[best]=0
  if(sum(ben.rem)>0){
  newbens=ben.rem+ben[best]
  newcosts=cost.rem+cost[best]
  ranks=rank(newcosts/(newbens*pext),ties.method="first")
  best=unique(c(best,which(ranks==1)))
  cost.rem2=0.5*cost.rem
  cost.rem2[which(cands==clearing.ind)]=cost[which(cands==clearing.ind)]
  ben.rem[best]=cost.rem2[best]=0
  for(cr in 1:nrow(nogos)){
  if(sum(cands[best[2]]%in%nogos[cr,])>0){ben.rem[cands%in%nogos[cr,]]=0}}
  newbens=ben.rem+sum(ben[best])
  if(sum(newbens)>0){
  if(clearing.ind%in%cands[best]){pext=1}
  newcosts=cost.rem2+cost[best[1]]+ifelse(length(best)==2,cost.rem[best[2]],0)
  ranks=rank(newcosts/(newbens*pext),ties.method="first")
  best=sort(unique(c(best,which(ranks==1))))
  }}}
  multiplier=c(1,0.5,0.25)[1:length(best)]
  best0=cands[best]
  bestcol=ifelse(length(best0)==1,best0,999)
  tben=sum((ben*pext)[best],na.rm=T)
  tcost=sum(cost[best]*multiplier,na.rm=T)
  bestid=paste(best0,collapse="+")
}
  return(list(bestcol=bestcol,bestid=bestid,ben=tben,cost=tcost))  
}

#argmax.comb=cmpfun(argmax.comb)


#get combos

namethem=function(cc,actnames){
  if(length(grep("\\+",cc))==0){vect=actnames[as.numeric(cc)]}else{
  vect=actnames[sort(as.numeric(strsplit(cc,"\\+")[[1]]))]}
  return(c(vect,rep(NA,3-length(vect))))
}

#calc species specific benefit for combo actions...

combben=function(bens,donothing){
  
}



argmax.comb.shiny=function(progind=1,ben,cost, clearing.prob=0,clearing.ind=1,all=999,nogos=nogo.combos){
  bestcol=bestid=tben=tcost=NA
  cands=which(ben/cost>0)
  pext=rep(1-clearing.prob,length(ben))
  pext[clearing.ind]=1
  if(length(cands)>0){
    ben=ben[cands]
    cost=cost[cands]
    pext=pext[cands]
    ranks=rank(cost/(ben*pext),ties.method="first")
    best=which(ranks==1)
    best.val=((ben*pext)/cost)[best]
    if(cands[best]==clearing.ind){pext=1}
    if(cands[best]!=all){
      if(all%in%cands){ben[cands==all]=0}
      for(cr in 1:nrow(nogos)){
        if(sum(cands[best]%in%nogos[cr,])>0){ben[cands%in%nogos[cr,]]=0}}
      ben.rem=ben
      cost.rem=0.5*cost
      cost.rem[which(cands==clearing.ind)]=cost[which(cands==clearing.ind)]
      ben.rem[best]=cost.rem[best]=0
      if(sum(ben.rem)>0){
        newbens=ben.rem+ben[best]
        newcosts=cost.rem+cost[best]
        ranks=rank(newcosts/(newbens*pext),ties.method="first")
        best=unique(c(best,which(ranks==1)))
        cost.rem2=0.5*cost.rem
        cost.rem2[which(cands==clearing.ind)]=cost[which(cands==clearing.ind)]
        ben.rem[best]=cost.rem2[best]=0
        for(cr in 1:nrow(nogos)){
          if(sum(cands[best[2]]%in%nogos[cr,])>0){ben.rem[cands%in%nogos[cr,]]=0}}
        newbens=ben.rem+sum(ben[best])
        if(sum(newbens)>0){
          if(clearing.ind%in%cands[best]){pext=1}
          newcosts=cost.rem2+cost[best[1]]+ifelse(length(best)==2,cost.rem[best[2]],0)
          ranks=rank(newcosts/(newbens*pext),ties.method="first")
          best=sort(unique(c(best,which(ranks==1))))
        }}}
    multiplier=c(1,0.5,0.25)[1:length(best)]
    best0=cands[best]
    bestcol=ifelse(length(best0)==1,best0,999)
    tben=sum((ben*pext)[best],na.rm=T)
    tcost=sum(cost[best]*multiplier,na.rm=T)
    bestid=paste(best0,collapse="+")
  }
  setProgress(value=progind)
  return(list(bestcol=bestcol,bestid=bestid,ben=tben,cost=tcost))  
}

##map specific actions (mainly for shiny)

showaction<-function(showact,bestact,lookup){
  showcodes<-lookup$actcodes[lookup$X1%in%showact | lookup$X2%in%showact | lookup$X3%in%showact]
  showras<-bestact%in%showcodes
  values(showras)=ifelse(values(showras)==1,1,NA)
  return(showras)
}


