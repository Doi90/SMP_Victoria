################################################
################################################
###                                          ###
### CONVERT PRE-1750 FLORA RASTERS TO EXTANT ###
###                                          ###
###   This script converts the pre-1750      ###
### flora species rasters into their extant  ###
### extent.                                  ###
###                                          ###
################################################
################################################

## Load packages

library(raster)
library(stringr)

## Load current native vegetation extent mask

nv_mask <- raster("data/masks/nv225vg.tif")

## Load sample raster with correct CRS

ras_CRS <- raster("data/masks/vicmask225vg.tif")

## Extract all raster filenames

ras_list <- list.files("data/TO_BE_SORTED/SMP_Pre1750_HDMS_MultiSppFlora/",
                       pattern = ".tif$")

ras_list <- str_replace(ras_list,
                        "MultiSppSubSpp50RF66All_",
                        "")

ras_list2 <- list.files("data/sdm_outputs/",
                        pattern = ".tif$")

## Keep filenames of only flora that don't have a single-species
## model output already. This folder includes fauna as well, and they
## are aready done

ras_list <- ras_list[-which(ras_list %in% ras_list2)]

## Loop through all files to mask and save to correct place

for(i in seq_len(length(ras_list))){
  
  filepath_old <- sprintf("data/TO_BE_SORTED/SMP_Pre1750_HDMS_MultiSppFlora/MultiSppSubSpp50RF66All_%s",
                          ras_list[i])
  
  tmp <- raster(filepath_old)
  
  tmp <- mask(tmp, mask = nv_mask)
  
  projection(tmp) <- projection(ras_CRS)
  
  filepath_new <- sprintf("data/sdm_outputs/%s",
                          ras_list[i])
  writeRaster(tmp,
              filepath_new,
              overwrite = TRUE)

}
