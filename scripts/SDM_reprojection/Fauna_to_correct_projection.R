###################################################
###################################################
###                                             ###
### CONVERT FAUNA RASTERS TO CORRECT PROJECTION ###
###                                             ###
###   This script puts all fauna SDM rasters    ###
### in the correct projection                   ###
###                                             ###
###################################################
###################################################

## Load packages

library(raster)

## Load sample raster with correct CRS

ras_CRS <- raster("data/masks/vicmask225vg.tif")

## Extract all raster filenames

ras_list <- list.files("data/TO_BE_SORTED/SMP_HDMS_Fauna100percentMaskExtant/",
                       pattern = ".tif$")

## Loop through all files to mask and save to correct place

for(i in seq_len(length(ras_list))){
  
  filepath_old <- sprintf("data/TO_BE_SORTED/SMP_HDMS_Fauna100percentMaskExtant/%s",
                          ras_list[i])
  
  tmp <- raster(filepath_old)
  
  projection(tmp) <- projection(ras_CRS)
  
  filepath_new <- sprintf("data/sdm_outputs/%s",
                          ras_list[i])
  writeRaster(tmp,
              filepath_new,
              overwrite = TRUE)
  
}
