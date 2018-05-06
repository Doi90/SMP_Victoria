##########################################################
##########################################################
###                                                    ###
###                 ACTION BENEFITS                    ###
###                                                    ###
###	This script models the SMP expert elicitation (EE) ###
### responses as functions of species traits and       ###
### spatial context (threats, veg condition, etc). It  ###
### then creates rasters of expected action benefits   ###
### for each species * action combination.             ###
###                                                    ###
##########################################################
##########################################################


#####################
### Load Packages ###
#####################

library(mboost)
library(reshape2)
library(doParallel)
library(raster)
library(tcltk)

##############################
### Set-up Directory Paths ###
##############################

## Directory paths for EE data & SDMs

fdir_in_EE <- "data/expert_elicitation"

fdir_in_threats <- "data/threat_maps"

fdir_in_actions <- "data/action_costs"

fdir_in_SDM <- "data/sdm_outputs"

fdir_in_traits <- "data/trait_models"

fdir_in_masks <- "data/masks"

## Directory path for benefit raster outputs NB: will have to change for each simulation

fdir_out<-"outputs/benefits"

#################
### Load Data ###
#################

load("data/expert_elicitation/response.df.Rdf")
load("data/expert_elicitation/EEsitedata.Rdf")

predictor.stack <- stack(list.files("data/action_costs",
                                    pattern = ".tif$",
                                    full.names = TRUE))

response.df$group <- gsub("Bats",
                          "Mammals",
                          response.df$group)

rsp0 <- response.df

grpnms <- c("Mammal", "Bird", "Amphibian", "Reptile", "Plant")

mstopvals <- c(420, 640, 315, 500, 580)

mstopvals.na <- c(200, 300, 50, 100, 100)

#####################
### RUN THE MODEL ###
#####################

## Looping over taxon groups

for(k in seq_len(5)){  
  
  group <- sprintf("%ss", grpnms[k])
  
  responses0 <- rsp0[rsp0$group == group, ]
  
  load(sprintf("data/trait_models/%sTraits_EEModels.Rdf",
               grpnms[k]))
  
  ## Assign the taxon specific trait data to 'traits'
  
  if(k == 1){
    
    traits <- mamm.traits
    
  } else if(k == 2){
    
    traits <- bird.traits
    
  } else if(k == 3){
     
    traits <- frog.traits
     
  } else if(k == 4){
    
    traits <- reptile.traits
    
  } else {
     
    traits <- plant.traits
  }
   
  ## Model Benefits...currently uses meanBen...i.e. experts best estimate
  
  responses <- responses0[responses0$action != "Do Nothing", ] 
  
  threats <- sitedata[responses$ref.scenario, ]
  
  ## Get scale, strata, climate change settings from scenario title
  
  scale <- as.numeric(gsub(".*_", "", gsub("ha_.*", "", responses$title)))
  
  scen <- gsub("_st.*", "", responses$ref.scenario)
  
  strata <- as.numeric(gsub("_sc.*", "", gsub(".*_st", "", responses$ref.scenario)))
  
  nocc <- rep(0, length(strata))
  
  nocc[grep("_noCC_", responses$ref.scenario)] <- 1
  
  traitclass <- unlist(lapply(traits, class))
  
  for(t in names(traitclass)[traitclass == "integer"]){
    
    traits[t] <- as.numeric(traits[ , t])
    
  } 
  
  trait <- traits[match(responses$TaxonCode,
                        traits$TaxonCode),
                  -c(1:3)] #get traits (removing ID variables)
   
  ## Model training data
  
  mdat <- data.frame(y = as.numeric(responses$meanBen),
                     action = responses$action,
                     species = factor(responses$TaxonCode),
                     expert = responses$expert,
                     prob0 = responses$refmean/10,
                     tempscale = responses$tempscale,
                     condition = responses$condition,
                     scale,
                     threats,
                     strata,
                     nocc,
                     trait)
  
  mdat <- mdat[is.finite(mdat$y), ]
  
  mdat$species <- factor(mdat$species,
                         levels=c(levels(mdat$species),
                                  "newspp"))
   
  fmla <- sprintf("y ~ btree(%s) + brandom(expert)", # + brandom(scen)
                  paste(names(mdat)[!names(mdat) %in% c("y","expert")],
                        collapse=" , "))
  
  modelBen <- mboost(as.formula(fmla),
                     data = mdat,
                     control = boost_control(mstop = mstopvals[k],
                                             trace = TRUE))
  
  ## Predict persistence probs conditional on No Action.
  
  responses <- responses0[responses0$action == "Do Nothing", ]
  
  threats <- sitedata[responses$ref.scenario, ]
  
  ## Get scale, strata, climate change settings from scenario title
  
  scale <- as.numeric(gsub(".*_", "", gsub("ha_.*", "", responses$title)))
  
  scen <- gsub("_st.*", "", responses$ref.scenario)
  
  strata <- as.numeric(gsub("_sc.*", "", gsub(".*_st", "", responses$ref.scenario)))
  
  nocc <- rep(0, length(strata))
  
  nocc[grep("_noCC_", responses$ref.scenario)] <- 1
  
  trait <- traits[match(responses$TaxonCode,
                        traits$TaxonCode),
                  -c(1:3)]
  
  ## Model training data
  
  mdat.na <- data.frame(y = as.numeric(responses$refmean),
                        action = responses$action,
                        species = factor(responses$TaxonCode),
                        expert = responses$expert,
                        tempscale = responses$tempscale,
                        condition = responses$condition,
                        scale,
                        threats,
                        strata,
                        nocc,
                        trait)
  
  mdat.na <- mdat.na[is.finite(mdat.na$y), ]
  
  mdat.na$species <- factor(mdat.na$species,
                            levels = c(levels(mdat.na$species),
                                       "newspp"))
  
  #mdat.na$expert=factor(mdat.na$expert,levels=c(unique(mdat.na$expert),"AAA"))
  
  ## Build model
  
  fmla.na <- sprintf("y ~ btree(%s) + brandom(expert)", # + brandom(scen)
                     paste(names(mdat.na)[!names(mdat.na) %in% c("y","expert")],
                           collapse=" , "))
  
  modelDN <- mboost(as.formula(fmla.na),
                    data = mdat.na,
                    control = boost_control(mstop = mstopvals.na[k],
                                            trace = TRUE))
  
  ##SAVE MODEL OBJECTS (if mapping predictions separately)
  
  #save(modelBen,file=file.path(fdir_out,"modelBen.Rmod"))
  #save(modelDN,file=file.path(fdir_out,"modelDN.Rmod"))
  
  
  ##MAP MODEL PREDICTIONS ACROSS EXTENT OF HDMS..... 
  
  
  savenames <- c("Harvesting",
                 "Weeds",
                 "Rabbits",
                 "Dom.Grazing",
                 "Cats",
                 "Foxes",
                 "Fox.Cats",
                 "Graz.Weeds",
                 "Reveg",
                 "Deer",
                 "Fuelred",
                 "Goat",
                 "Horses",
                 "Miners",
                 "Phytoph",
                 "Pigs")
  
  actionnames <- c("StopHarvesting",
                   "WeedControl",
                   "RabbitControl",
                   "StockExclusion",
                   "CatControl",
                   "FoxControl",
                   "FoxandCatControl",
                   "WeedControlandStockExclusion",
                   "Reveg",
                   "DeerControl",
                   "ReducedFRBfrequency",
                   "FeralGoatControl",
                   "FeralHorseControl",
                   "MinerControl",
                   "Phytophthoramanagement",
                   "FeralPigControl")
  
  vic <- raster("data/masks/vicmask225vg.tif")
    
  vicnv <- raster("data/masks/nv225vg.tif")
  
  threat.action <- list(DeerControl = "deer",
                        FeralGoatControl = "goat",
                        FeralPigControl = "pigs",
                        Phytophthoramanagement = "phytophthora",
                        RabbitControl = "rabbit",
                        ReducedFRBfrequency = "fuelburning",
                        StockExclusion = "grazing",
                        StopHarvesting = "harvesting",
                        WeedControl = "transweed",
                        WeedControlandStockExclusion = c("grazing","transweed"))
  
  traitnames <- names(mdat)[names(mdat) %in% names(traits)]
  
  tscale <- 50
  
  scale <- 100
  
  nocc <- 0
   
  nvcells <- which(values(vicnv) == 1)
  
  predictor.vals <- values(predictor.stack)
  
  # if(group != "Plants"){
  #   
  #   hdmlist <- list.files("O:/SAN_Projects/NPRINT4/HDMV52016/Tiff225/Masked/SMP_HDMS_Fauna100percentMaskExtant/",pattern=".tif$",full=T)
  # }else{
  #   hdmlist=list.files("O:/SAN_Projects/NPRINT4/HDMV52016/Tiff225/Masked/SMP_HDMS/",pattern=".tif$",full=T)
  #   hdmlist2=list.files("O:/SAN_Projects/NPRINT4/HDMV52016/Tiff225/Masked/SMP_HDMS_MultiSppFlora",pattern=".tif$",full=T)
  #   hdmlist=c(hdmlist,hdmlist2)  
  # }
  
  hdmlist <- list.files("data/sdm_outputs/",
                        pattern = ".tif$")
  
  mapacts <- levels(factor(responses0$action))
  
  mapacts <- c("Do Nothing", 
               mapacts[mapacts != "Do Nothing"])
  
  mapacts <- mapacts[!mapacts %in% c("Phytophthora management",
                                     "Reduced FRB frequency",
                                     "Stop Harvesting")]
   
  ## Imput missing trait values (assuming model used can't handle NA's
  ## for prediction)
  
  if(!"TaxonGroup" %in% names(traits)){
    
    traits$TaxonGroup <- 1
    
  } 
  
  traits.impt <- traits
  
  meanfun <- function(x){
    
    if(is.numeric(x)){
      
      x <- ifelse(is.na(x),
                  mean(x, na.rm = TRUE),
                  x)
      
    } else { 
      
      levs <- levels(x)
      
      meanval <- names(sort(table(x),
                            decreasing = TRUE))[1]
      
      x[is.na(x)] <- meanval
      
    } 
    
    return(x)
    
  }
   
  if(any(is.na(traits$TaxonGroup))){
    
    f2 <- factor(traits$TaxonGroup,
                 levels = c(levels(traits$TaxonGroup),
                            "MISSING"))
    
    f2[is.na(f2)] <- "MISSING"
    
    traits$TaxonGroup <- f2
    
  } 
  
  avgover <- traits$TaxonGroup
  
  if(group == "Plants"){
    
    avgover = traits$florafauna_spp_lu_Family
    
  } 
  
  ## Impute trait values for species with missing values
  
  for(c in 4:ncol(traits)){
    
    traits.impt[ , c] <- unlist(tapply(traits.impt[ , c],
                                       avgover,
                                       meanfun))
  } 
  
  NoRespSpp <- traits$TaxonCode[traits$TaxonCode %in% responses$TaxonCode == FALSE]
  
  doSpecies <- traits$TaxonCode[1:5]
  
  ## Set up clusters for running in parallel
  
  nnodes <- detectCores(logical = TRUE) - 1
  
  cl <- makeCluster(nnodes)
  
  registerDoParallel(cl)
  
  foreach(sp = iter(doSpecies),
          # .packages = c("raster","tcltk"),
          .packages = c("raster"),
          .errorhandling = "stop") %dopar% {
             
            usehdm <- grep(sprintf("Spp%s", sp),
                           hdmlist)
            
            if(length(usehdm) > 0){
              
              ## Load raster, covert back to probability
              ## Pretty sure the [1] is not required, but doesn't hurt to leave in
              
              hdm <- raster(sprintf("data/sdm_outputs/%s",
                                    hdmlist[usehdm[1]])) / 100
              
              habvals <- which(values(hdm) > 0)
              
              # pbsp <- tkProgressBar(paste("Mapping benefits for", sp),
              #                       min = 0,
              #                       max = nlevels(mdat$action),
              #                       label = paste(length(mapacts),"actions to do"))
              
              for(actn in mapacts){
                
                if(sp %in% NoRespSpp){
                  
                  possrange <- quantile(mdat$y,
                                        c(0.1, .9))
                  
                } else {
                  
                  possrange <- range(mdat$y[mdat$species == sp & mdat$action == actn],
                                     na.rm = TRUE) / 10
                  
                  if(sum(is.finite(possrange)) < 2){
                    
                    possrange <- range(mdat$y[mdat$species == sp], 
                                       na.rm = TRUE) / 10
                    
                  }
                  
                  if(sum(is.finite(possrange)) < 2){
                    
                    possrange <- quantile(mdat$y,
                                          c(0.1,.9))
                    
                  }
                }
                
                usemod <- modelBen
                usemdat <- mdat
                
                if(actn == "Do Nothing"){
                  
                  usemod <- modelDN
                  usemdat <- mdat.na
                  
                }
                
                usevals <- habvals
                aname <- gsub(" ",
                              "",
                              actn)
                
                if(aname %in% names(threat.action)){
                  
                  tvals <- data.frame(predictor.vals[habvals,threat.action[[aname]]])
                  
                  usevals <- habvals[which(apply(tvals, 1, prod) > 0)]
                  
                }
                
                nvals <- length(usevals)
                
                preds <- vic * 0
                
                #values(preds)=ifelse(values(vic)==1,0,NA)
                
                if(nvals > 10){
                  
                  const <- data.frame(species = factor(sp,
                                                       levels = levels(usemdat$species)),
                                      action = factor(actn,
                                                      levels = levels(usemdat$action)),
                                      tempscale = tscale,
                                      scale = scale,
                                      nocc = 0)
                  
                  const <- data.frame(const,
                                      traits.impt[match(sp,
                                                        traits.impt$TaxonCode),
                                                  traitnames])
                  
                  if(sp %in% NoRespSpp){
                    
                    const$species <- factor("newspp",
                                            levels = levels(usemdat$species))
                    
                  }
                  
                  newdat <- data.frame(predictor.vals[usevals, ],
                                       const)
                  
                  if(actn != "Do Nothing"){
                    
                    newdat$prob0 <- donothing[usevals]
                    
                  }
                  
                  ##THIS SPLITS RASTER UP INTO 200K cell chunks for processing...loops over chunks..PROB NOT EFFICIENT!
                  
                  nstarts <- ceiling(nvals / 200000)
                  
                  starts <- seq(1, 
                                by = 200000,
                                length.out = nstarts)
                  
                  ends <- c(starts[-1] - 1,
                            nvals)
                  
                  # pb <- tkProgressBar(paste("Predicting",
                  #                           actn),
                  #                     min = 0,
                  #                     max = nstarts,
                  #                     label = paste(nstarts,
                  #                                   "subsets to do"))
                  
                  for(ind in seq_len(nstarts)){
                    
                    newvals <- predict(usemod,
                                       newdata = newdat[starts[ind]:ends[ind], ],
                                       which = 1) + usemod$offset
                    
                    values(preds)[usevals[starts[ind]:ends[ind]]] <- newvals / 10
                    
                    # setTkProgressBar(pb,
                    #                  ind,
                    #                  label = paste(ind,
                    #                                "of",
                    #                                nstarts,
                    #                                "subsets done"))
                    
                    gc()
                    
                  }
                  
                  #close(pb)
                  
                  if(actn == "Do Nothing"){
                    
                    values(preds) <- ifelse(values(preds) < 0,
                                            0,
                                            ifelse(values(preds) > 1,
                                                   1,
                                                   values(preds)))
                    
                    donothing <- values(preds)
                    
                    preds <- preds * hdm
                    
                    DNpreds <- values(preds)
                    
                  }
                  
                  if(actn != "Do Nothing"){
                    
                    preds <- preds * hdm
                    
                    values(preds) <- ifelse(values(preds) > (1 - DNpreds),
                                            (1 - DNpreds),
                                            ifelse(values(preds) < (-DNpreds),
                                                   -DNpreds,
                                                   values(preds)))
                    
                  }
                  
                }
                
                if(actn != "Do Nothing"){
                
                shortname <- savenames[match(aname,
                                             actionnames)]
                
                tifname <- sprintf("Spp%s_delta_%s.tif",
                                   sp,
                                   shortname) 
                }
                
                if(actn == "Do Nothing"){
                  
                  tifname <- sprintf("Spp%2_DoNothing.tif",
                                     sp)
                }
                
                writeRaster(preds,
                            file = sprintf("outputs/benefits/%s",
                                           tifname),
                            format = "GTiff",
                            overwrite = TRUE)
                
                # setTkProgressBar(pbsp,
                #                  match(actn,
                #                        mapacts),
                #                  label = paste(actn,
                #                                "done")) 
                
                
              }
              
              # close(pbsp)
              
              gc()
              
            }
            
            # if(!exists("pb0")){
            #   
            #   pb0 <- tkProgressBar("Mapping Species Benefits",
            #                        min = 0,
            #                        max = length(doSpecies))
            #   
            #   setTkProgressBar(pb0,
            #                    match(sp,
            #                          doSpecies))
            #   
            # }
            
            stopCluster(cl)
            
            gc()
            
            
          }
}














