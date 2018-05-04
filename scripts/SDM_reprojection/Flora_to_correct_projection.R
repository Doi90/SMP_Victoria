###################################################
###################################################
###                                             ###
### CONVERT FLORA RASTERS TO CORRECT PROJECTION ###
###                                             ###
###   This script puts all flora SDM rasters    ###
### in the correct projection                   ###
###                                             ###
###################################################
###################################################

## Load packages

library(raster)

## Load sample raster with correct CRS

ras_CRS <- raster("data/masks/vicmask225vg.tif")

## Extract all raster filenames

ras_list <- list.files("data/TO_BE_SORTED/SMP_HDMS/",
                       pattern = ".tif$")

ras_list2 <- list.files("data/TO_BE_SORTED/SMP_HDMS_Fauna100percentMaskExtant",
                       pattern = ".tif$")

## Keep filenames of only flora
## This folder includes fauna as well, and they
## are aready done

ras_list <- ras_list[-which(ras_list %in% ras_list2)]

## Loop through all files to mask and save to correct place

for(i in seq_len(length(ras_list))){
  
  filepath_old <- sprintf("data/TO_BE_SORTED/SMP_HDMS/%s",
                          ras_list[i])
  
  tmp <- raster(filepath_old)
  
  projection(tmp) <- projection(ras_CRS)
  
  filepath_new <- sprintf("data/sdm_outputs/%s",
                          ras_list[i])
  writeRaster(tmp,
              filepath_new,
              overwrite = TRUE)
  
}
