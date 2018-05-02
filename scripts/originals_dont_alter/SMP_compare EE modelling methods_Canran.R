#	This program is to do sensitivity analysis for SMP.

#	In this program we model the two probabilities (prob_under_action and prob_ref) separately first,
#	then take their difference as the benifit. 
#	We also model the benifit at the same time.
#	In this way, it's better for comparing the two approaches.

#	In this program we also use REEMtree (can only take one random effect variable: species because it can only take one random effect) 
#	in addition to random forest and mboost. 

#	This program is similar to SMP-7.r. But in SMP-7.r, the correspondence between reference and action
#	seems have a problem. This was corrected here.  

#	Another difference is that the trait florafauna_spp_lu_Family for plants and TaxonGrp for birds
#	are not used for RF models since they have too many categories.

#	Use load() to load the R objects.

#	source("\\\\ari-vm-00\\SpatialStorage\\SAN_Projects\\liu\\U_copy\\accudx\\SMP-8.r")



set.seed(677007)


#----------------------------------------------------------------------------------------------------

#	calculate accuracy measures

accudx <- function(y0,y1) {
	bias <- mean(y1-y0)
	abs_bias <- mean(abs(y0-y1))
	rmse <- sqrt(mean((y0-y1)^2))

	c(bias,abs_bias,rmse)
   
}


#----------------------------------------------------------------------------------------------------

#	calculating aberrance

library(REEMtree)
library(randomForest)
library(mboost)
library(openxlsx)
library(reshape2)
#library(doParallel)

rk <- 7		# just to differentiate different runs of results

nsm <- 100

na0=function(x){ifelse(is.na(as.numeric(x)),0,as.numeric(x))}
sumna0=function(x){sum(na0(x))}



fdir_in <- "\\\\ari-vm-00\\SpatialStorage\\SAN_Projects\\liu\\U_copy\\stat_consult\\smp\\data\\"
fdir_out <- "\\\\ari-vm-00\\SpatialStorage\\SAN_Projects\\liu\\U_copy\\stat_consult\\smp\\rst8\\"

#fdir_in <- "S:/SAN_Projects/liu/U_copy/stat_consult/smp/data/"
#fdir_out <- "S:/SAN_Projects/liu/U_copy/stat_consult/smp/rst8/"



load(paste(fdir_in,"response.df.Rdf",sep=""))
rsp0  = response.df

fln_site <- paste(fdir_in,"EEresults_all_13October2015.xlsx",sep="")
sitedata=read.xlsx(fln_site,sheet=2,rowNames=T)
sitedata[,1:ncol(sitedata)]=apply(sitedata[,1:ncol(sitedata)],2,na0)

grpnms = c("Mammal","Plant","Bird","Amphibian","Reptile")
trtnms = c("Mammal","Plant","Bird","Frog","Reptile")


for (k in c(5)) {

group = paste(grpnms[k],"s",sep="")

responses0 = rsp0[rsp0$group==group,c(1:12,15,18,20)]

load(paste(fdir_in,trtnms[k],"Traits_EEModels.R",sep=""))

if (k==1) {
	traits = mamm.traits

##### here I add a few lines to deal with the messy data in traits$nesting.code because there are both "N" and "N ", they should be the same, also change "" to NA (it is missing data). 

levels(traits$Diet)[14] <- "O"		# it originally was "O "
levels(traits$Diet)[8] <- "I+N"		# it originally was "I +N"
levels(traits$Diet)[7] <- "I"		# it originally was "I "
levels(traits$Diet)[1] <- NA		# it originally was ""

levels(traits$nesting.code)[10] <- "N"	# it originally was "N "
levels(traits$nesting.code)[1] <- NA	# it originally was ""

##### end


} else {
	if (k==2) {
		traits = traits
	} else {
		if (k==3) {
			traits = bird.traits
		} else {
			if (k==4) {
				traits = frog.traits
				names(traits)[1] <- "TaxonCode"
			} else {
				traits = reptile.traits
				names(traits)[1] <- "TaxonCode"
				levels(traits)[1] <- NA	# it originally was ""
			}
		}
	}
}



responses = responses0[responses0$action != "Do Nothing",]

threats=sitedata[responses$ref.scenario,]

#get scale, strata, climate change settings from scenario title

scale=as.numeric(gsub(".*_","",gsub("ha_.*","",responses$title)))
scen=gsub("_st.*","",responses$ref.scenario)
strata=as.numeric(gsub("_sc.*","",gsub(".*_st","",responses$ref.scenario)))
nocc=rep(0,length(strata))
nocc[grep("_noCC_",responses$ref.scenario)]=1

trait=traits[match(responses$TaxonCode,traits$TaxonCode),-c(1:3)]

mdat=data.frame(y=as.numeric(responses$mean),action=responses$action,tempscale=responses$tempscale,scale,threats,strata,nocc,trait)

mdat0 = mdat

fm_tmp = paste("y~",paste(names(mdat)[2:ncol(mdat)],collapse="+"))

fmla_mf=as.formula(fm_tmp)	# for REEMtree model

if (k==2) fm_tmp = paste("y~",paste(setdiff(names(mdat)[2:ncol(mdat)],"florafauna_spp_lu_Family"),collapse="+"))	# remove florafauna_spp_lu_Family from the formula since it has more than 100 categories
if (k==3) fm_tmp = paste("y~",paste(setdiff(names(mdat)[2:ncol(mdat)],"TaxonGrp"),collapse="+"))

fmla_rf=as.formula(fm_tmp)	# for random forest model


fmla=paste("y~btree(",paste(names(mdat)[2:ncol(mdat)],collapse=","),")")		# for boost model
fmla2=gsub("action","species, action",fmla)


fmla=as.formula(paste(fmla,"+ brandom(species) + brandom(expert)"))  #+ brandom(scen) 
fmla2=as.formula(paste(fmla2," + brandom(expert)")) # + brandom(scen)

mdat=data.frame(mdat,expert=responses$expert,scen,species=factor(responses$TaxonCode),ref.scenario=responses$ref.scenario)



#First model (fmla) has random species effects...used to predict to species not included in oringal expert elicitation (based on traits only) 

mdat <- na.omit(mdat)

mdat = unique(mdat)

mod=mboost(fmla,data=mdat,control=boost_control(mstop=1000,trace=T))

#cvm <- cvrisk(mod)

#mod[mstop(cvm)]


#get fitted values
preds=predict(mod,newdata=mdat,which=1)
preds=preds+attr(preds,"offset")
plot(preds~factor(mdat$y))

dx0 <- accudx(preds,mdat$y) 


#--------------------------------------------------------------------------------------------------

###NOW directly predict benifit.

#	here, fmla_rf is the same as in the above.

mdat.bft=mdat0

mdat.bft$y = as.numeric(responses$meanBen)

mdat.bft=data.frame(mdat.bft,expert=responses$expert,scen,species=factor(responses$TaxonCode),ref.scenario=responses$ref.scenario)

# random species model

mdat.bft <- na.omit(mdat.bft)
mdat.bft = unique(mdat.bft)

mod.bft=mboost(fmla,data=mdat.bft,control=boost_control(mstop=1000,trace=T))

#plot(fitted(mod.bft)~factor(mdat.bft$y))	

#cvm <- cvrisk(mod.bft)

#mod.na[mstop(cvm)]

preds=predict(mod.bft,newdata=mdat.bft,which=1)
preds=preds+attr(preds,"offset")
plot(preds~factor(mdat.bft$y))

dx0 <- accudx(preds,mdat.bft$y) 


#--------------------------------------------------------------------------------------------------

###NOW predict persistence probs conditional on No Action.

responses = responses0[responses0$action == "Do Nothing",]

threats=sitedata[responses$ref.scenario,]

#get scale, strata, climate change settings from scenario title

scale=as.numeric(gsub(".*_","",gsub("ha_.*","",responses$title)))
scen=gsub("_st.*","",responses$ref.scenario)
strata=as.numeric(gsub("_sc.*","",gsub(".*_st","",responses$ref.scenario)))
nocc=rep(0,length(strata))
nocc[grep("_noCC_",responses$ref.scenario)]=1

trait=traits[match(responses$TaxonCode,traits$TaxonCode),-c(1:3)]

mdat.na = data.frame(y=as.numeric(responses$mean),tempscale=responses$tempscale,scale,threats,strata,nocc,trait)

mdat.na$y = as.numeric(responses$refmean)


fm_tmp = paste("y~",paste(names(mdat.na)[2:ncol(mdat.na)],collapse="+"))

fmla.na_mf=as.formula(fm_tmp)	# for REEMtree model


if (k==2) fm_tmp = paste("y~",paste(setdiff(names(mdat.na)[2:ncol(mdat.na)],"florafauna_spp_lu_Family"),collapse="+"))	# remove florafauna_spp_lu_Family from the formula since it has more than 100 categories
if (k==3) fm_tmp = paste("y~",paste(setdiff(names(mdat.na)[2:ncol(mdat.na)],"TaxonGrp"),collapse="+"))

fmla.na_rf=as.formula(fm_tmp)	# for random forest model

fmla.na=paste("y~btree(",paste(names(mdat.na)[2:ncol(mdat.na)],collapse=","),")")
fmla2.na=gsub("tempscale","species, tempscale",fmla.na)
fmla.na=as.formula(paste(fmla.na,"+ brandom(species) + brandom(expert)"))  #+ brandom(scen) 
fmla2.na=as.formula(paste(fmla2.na," + brandom(expert)")) # + brandom(scen)

mdat.na=data.frame(mdat.na,expert=responses$expert,scen,species=factor(responses$TaxonCode),ref.scenario=responses$ref.scenario)



# random species model

mdat.na <- na.omit(mdat.na)
mdat.na = unique(mdat.na)

mod.na=mboost(fmla.na,data=mdat.na,control=boost_control(mstop=1000,trace=T))

#plot(fitted(mod.na)~factor(mdat.na$y))	

#cvm <- cvrisk(mod.na)

#mod.na[mstop(cvm)]

preds=predict(mod.na,newdata=mdat.na,which=1)
preds=preds+attr(preds,"offset")
plot(preds~factor(mdat.na$y))

dx0 <- accudx(preds,mdat.na$y) 




#------------------------------------------------------------------------------

nt1 <- dim(mdat)[1]
nt2 <- dim(mdat.na)[1]
nt4 <- dim(mdat.bft)[1]

all_spp <- intersect(unique(mdat$species),unique(mdat.na$species))
all_spp <- intersect(all_spp,unique(mdat.bft$species))

nspp <- length(all_spp)

m <- floor(nspp/4)


dx1 <- NULL	# for accuracy index
dx2 <- NULL
dx3 <- NULL
dx4 <- NULL

dt3 <- NULL
dt4 <- NULL


for (i in 1:nsm) {

cat("sm ",i," start","\n")

#	for prob of persistence under action

	tst_spp <- sample(all_spp,m)

	id <- mdat$species %in% tst_spp

	n1 <- sum(id)

	mdat_tst <- mdat[id,]

	mdat_trn <- mdat[!id,]

#	mboost model

	mod=mboost(fmla,data=mdat_trn,control=boost_control(mstop=1000,trace=T))

	if (k != 2) {

	cvm <- cvrisk(mod)

	mod[mstop(cvm)]

	}

	preds=predict(mod,newdata=mdat,which=1)
	preds=preds+attr(preds,"offset")



#	random forest model

	rf <- randomForest(fmla_rf, data=mdat_trn)
	pred = predict(rf,mdat)



#	merf model


	mf <-REEMtree(fmla_mf,data=mdat_trn,random=~1|species)
	prd = predict(mf, mdat, id=mdat$species,EstimateRandomEffects=TRUE)



	dx1 <- rbind(dx1,c(i,1,0,accudx(preds[!id],mdat$y[!id])))	# 1: mboost, 0: training data 
	dx1 <- rbind(dx1,c(i,1,1,accudx(preds[id],mdat$y[id])))		# 1: mboost, 1: test data  
	dx1 <- rbind(dx1,c(i,2,0,accudx(pred[!id],mdat$y[!id])))	# 2: random forest, 0: training data 
	dx1 <- rbind(dx1,c(i,2,1,accudx(pred[id],mdat$y[id])))		# 2: random forest, 1: test data
	dx1 <- rbind(dx1,c(i,3,0,accudx(prd[!id],mdat$y[!id])))		# 3: merf, 0: training data 
	dx1 <- rbind(dx1,c(i,3,1,accudx(prd[id],mdat$y[id])))		# 3: merf, 1: test data


	tmp11 <- data.frame(sm=rep(i,nt1-n1),trn_tst=rep(0,nt1-n1),species=mdat$species[!id],ref.scenario=mdat$ref.scenario[!id],
			expert=mdat$expert[!id],prob_obs=mdat$y[!id],prob_pred_mb=preds[!id],
			prob_pred_rf=pred[!id],prob_pred_mf=prd[!id])
	tmp12 <- data.frame(sm=rep(1,n1),trn_tst=rep(1,n1),species=mdat$species[id],ref.scenario=mdat$ref.scenario[id],
			expert=mdat$expert[id],prob_obs=mdat$y[id],prob_pred_mb=preds[id],
			prob_pred_rf=pred[id],prob_pred_mf=prd[id])


###	for prob of persistence with no action

	id <- mdat.na$species %in% tst_spp

	n2 <- sum(id)

	mdat.na_tst <- mdat.na[id,]

	mdat.na_trn <- mdat.na[!id,]


#	mboost model

	mod=mboost(fmla.na,data=mdat.na_trn,control=boost_control(mstop=1000,trace=T))

	if (k != 2) {

	cvm <- cvrisk(mod)

	mod[mstop(cvm)]

	}


	preds=predict(mod,newdata=mdat.na,which=1)
	preds=preds+attr(preds,"offset")


#	random forest model

	rf <- randomForest(fmla.na_rf, data=mdat.na_trn)
	pred = predict(rf,mdat.na)


#	merf model


	mf <-REEMtree(fmla.na_mf,data=mdat.na_trn,random=~1|species)
	prd = predict(mf, mdat.na, id=mdat.na$species,EstimateRandomEffects=TRUE)



	dx2 <- rbind(dx2,c(i,1,0,accudx(preds[!id],mdat.na$y[!id])))	# 1: mboost, 0: training data 
	dx2 <- rbind(dx2,c(i,1,1,accudx(preds[id],mdat.na$y[id])))	# 1: mboost, 1: test data  
	dx2 <- rbind(dx2,c(i,2,0,accudx(pred[!id],mdat.na$y[!id])))	# 2: random forest, 0: training data 
	dx2 <- rbind(dx2,c(i,2,1,accudx(pred[id],mdat.na$y[id])))	# 2: random forest, 1: test data
	dx2 <- rbind(dx2,c(i,3,0,accudx(prd[!id],mdat.na$y[!id])))		# 3: merf, 0: training data 
	dx2 <- rbind(dx2,c(i,3,1,accudx(prd[id],mdat.na$y[id])))		# 3: merf, 1: test data
 

	tmp21 <- data.frame(sm=rep(i,nt2-n2),trn_tst=rep(0,nt2-n2),species=mdat.na$species[!id],ref.scenario=mdat.na$ref.scenario[!id],
				expert=mdat.na$expert[!id],prob_obs_ref=mdat.na$y[!id],prob_pred_mb_ref=preds[!id],
				prob_pred_rf_ref=pred[!id],prob_pred_mf_ref=prd[!id])
	tmp22 <- data.frame(sm=rep(i,n2),trn_tst=rep(1,n2),species=mdat.na$species[id],ref.scenario=mdat.na$ref.scenario[id],
				expert=mdat.na$expert[id],prob_obs_ref=mdat.na$y[id],prob_pred_mb_ref=preds[id],
				prob_pred_rf_ref=pred[id],prob_pred_mf_ref=prd[id])



#	here calculate benifit from the above two predictions 

	tmp1 <- merge(tmp11,tmp21,by=c("trn_tst","species","expert","ref.scenario"))

	tmp1$sm <- rep(i,dim(tmp1)[1])

	tmp2 <- merge(tmp12,tmp22,by=c("trn_tst","species","expert","ref.scenario"))

	tmp2$sm <- rep(i,dim(tmp2)[1])

	

	dx3 <- rbind(dx3,c(i,1,0,accudx(tmp1$prob_pred_mb-tmp1$prob_pred_mb_ref,tmp1$prob_obs-tmp1$prob_obs_ref)))
	dx3 <- rbind(dx3,c(i,1,1,accudx(tmp2$prob_pred_mb-tmp2$prob_pred_mb_ref,tmp2$prob_obs-tmp2$prob_obs_ref)))
	dx3 <- rbind(dx3,c(i,2,0,accudx(tmp1$prob_pred_rf-tmp1$prob_pred_rf_ref,tmp1$prob_obs-tmp1$prob_obs_ref)))
	dx3 <- rbind(dx3,c(i,2,1,accudx(tmp2$prob_pred_rf-tmp2$prob_pred_rf_ref,tmp2$prob_obs-tmp2$prob_obs_ref)))
	dx3 <- rbind(dx3,c(i,3,0,accudx(tmp1$prob_pred_mf-tmp1$prob_pred_mf_ref,tmp1$prob_obs-tmp1$prob_obs_ref)))
	dx3 <- rbind(dx3,c(i,3,1,accudx(tmp2$prob_pred_mf-tmp2$prob_pred_mf_ref,tmp2$prob_obs-tmp2$prob_obs_ref)))


	
	dt3 <- rbind(dt3,tmp1)
	dt3 <- rbind(dt3,tmp2)



#	directly model benifit

	id <- mdat.bft$species %in% tst_spp

	n4 <- sum(id)

	mdat.bft_tst <- mdat.bft[id,]

	mdat.bft_trn <- mdat.bft[!id,]


#	mboost model

	mod=mboost(fmla,data=mdat.bft_trn,control=boost_control(mstop=1000,trace=T))

	if (k != 2) {

	cvm <- cvrisk(mod)

	mod[mstop(cvm)]

	}

	preds=predict(mod,newdata=mdat.bft,which=1)
	preds=preds+attr(preds,"offset")


#	random forest model

	rf <- randomForest(fmla_rf, data=mdat.bft_trn)
	pred = predict(rf,mdat.bft)


#	merf model


	mf <-REEMtree(fmla_mf,data=mdat.bft_trn,random=~1|species)
	prd = predict(mf, mdat.bft, id=mdat.bft$species,EstimateRandomEffects=TRUE)


	dx4 <- rbind(dx4,c(i,1,0,accudx(preds[!id],mdat.bft$y[!id])))	# 1: mboost, 0: training data 
	dx4 <- rbind(dx4,c(i,1,1,accudx(preds[id],mdat.bft$y[id])))	# 1: mboost, 1: test data  
	dx4 <- rbind(dx4,c(i,2,0,accudx(pred[!id],mdat.bft$y[!id])))	# 2: random forest, 0: training data 
	dx4 <- rbind(dx4,c(i,2,1,accudx(pred[id],mdat.bft$y[id])))	# 2: random forest, 1: test data
	dx4 <- rbind(dx4,c(i,3,0,accudx(prd[!id],mdat.bft$y[!id])))	# 3: merf, 0: training data 
	dx4 <- rbind(dx4,c(i,3,1,accudx(prd[id],mdat.bft$y[id])))	# 3: merf, 1: test data
 


	dt4 <- rbind(dt4,data.frame(sm=rep(i,nt4-n4),trn_tst=rep(0,nt4-n4),species=mdat.bft$species[!id],ref.scenario=mdat.bft$ref.scenario[!id],
				expert=mdat.bft$expert[!id],bft_obs=mdat.bft$y[!id],bft_pred_mb=preds[!id],
				bft_pred_rf=pred[!id],bft_pred_mf=prd[!id]))
	dt4 <- rbind(dt4,data.frame(sm=rep(i,n4),trn_tst=rep(1,n4),species=mdat.bft$species[id],ref.scenario=mdat.bft$ref.scenario[id],
				expert=mdat.bft$expert[id],bft_obs=mdat.bft$y[id],bft_pred_mb=preds[id],
				bft_pred_rf=pred[id],bft_pred_mf=prd[id]))

	cat("sm ",i," finish","\n")
 
}


dx1 <- as.data.frame(dx1)
names(dx1) <- c("simulation","technique","data_trn0_tst1","mean_bias","mean_abs_bias","RMSE")

dx2 <- as.data.frame(dx2)
names(dx2) <- c("simulation","technique","data_trn0_tst1","mean_bias","mean_abs_bias","RMSE")

dx3 <- as.data.frame(dx3)
names(dx3) <- c("simulation","technique","data_trn0_tst1","mean_bias","mean_abs_bias","RMSE")

dx4 <- as.data.frame(dx4)
names(dx4) <- c("simulation","technique","data_trn0_tst1","mean_bias","mean_abs_bias","RMSE")


write.csv(dx1,file=paste(fdir_out,"predictive_accuracy_for_prob_under_action_",grpnms[k],"-",rk,".csv",sep=""),row.names=FALSE)
write.csv(dx2,file=paste(fdir_out,"predictive_accuracy_for_prob_ref_",grpnms[k],"-",rk,".csv",sep=""),row.names=FALSE)
write.csv(dx3,file=paste(fdir_out,"predictive_accuracy_benefit_indirect_",grpnms[k],"-",rk,".csv",sep=""),row.names=FALSE)
write.csv(dx4,file=paste(fdir_out,"predictive_accuracy_benefit_direct_",grpnms[k],"-",rk,".csv",sep=""),row.names=FALSE)

write.csv(dt3,file=paste(fdir_out,"predictions_persistence_prob_",grpnms[k],"-",rk,".csv",sep=""),row.names=FALSE)
write.csv(dt4,file=paste(fdir_out,"predictions_benefit_direct_",grpnms[k],"-",rk,".csv",sep=""),row.names=FALSE)




}





















