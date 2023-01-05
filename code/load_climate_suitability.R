## Code for loading koala suitability indices from Graham et al. 2019 for NSW properties into a CSV
#  Author: Frankie Cho
#  Created: 13 October 2022

library(tidyverse)
library(raster)
library(R.utils)
library(sf)
library(future)
library(future.apply)
library(exactextractr)
library(spatialEco)

# Modify "outpath" to choose where to save the output file
outpath <- "C:/Users/uqfcho/Documents/Koala/temporal_planning/"
#outpath <- '/Users/frankiecho/Documents/Github/koala-uncertainty/data/'

source('code/load_paths.R')

gcm <- c('cccma-cgcm31', 'ccsr-miroc32hi', 'ccsr-miroc32med', 'cnrm-cm3', 'csiro-mk30', 'gfdl-cm20', 
         'gfdl-cm21', 'giss-modeleh', 'giss-modeler', 'iap-fgoals10g', 'inm-cm30', 'ipsl-cm4',
         'mpi-echam5', 'mri-cgcm232a', 'ncar-ccsm30', 'ncar-pcm1', 'ukmo-hadcm3', 'ukmo-hadgem1')
rcp <- c('RCP45', 'RCP85')
years <- c('2015', '2025', '2035', '2045', '2055', '2065', '2075', '2085')

clim_scen <- expand.grid(rcp, gcm, years)
clim_scen_list <- split(clim_scen, seq(nrow(clim_scen))) %>%
  lapply(function(e) {
    names(e) <- c('rcp', 'gcm', 'year')
    e$id <- paste(e$rcp, e$gcm, e$year, sep='_')
    e$path <- paste(e$rcp, e$gcm, e$year, 'suitability', sep = '_')
    e$path <- paste0(rdm_path, suitability_path, e$path)
    if (file.exists(paste0(e$path, '.asc.gz'))) {
      e$path <- paste0(e$path, '.asc.gz')
      e$gz <- T
    } else {
      e$path <- paste0(e$path, '.asc')
      e$gz <- F
    }
    e
  })
names(clim_scen_list) <- lapply(clim_scen_list, function(x) x$id)

# Simple extraction based on point estimates
properties_source <- st_read(paste0(rdm_path, properties_path), layer = "spatial_pred_inperp_v1_0")
khab_thres <- raster(paste0(rdm_path, khab_path))
properties_khab <- exactextractr::exact_extract(khab_thres, properties_source, fun='sum')

properties <- properties_source %>%
  st_transform(4326)
properties_centroid <- properties_source %>%
  st_centroid() %>%
  st_transform(4326) %>%
  st_coordinates()

fcn_prop_suitability <- function(clim_scen_path, gz = T, polygon = T) {
  if (gz) {
    spdf <- gunzip(filename = clim_scen_path, remove = F)
  } else {
    spdf <- clim_scen_path
  }
  suit_raster <- spdf %>% 
    raster()
  if (polygon) {
    prop_suit <- raster::extract(suit_raster, properties)
  } else {
    prop_suit <- raster::extract(suit_raster, properties_centroid)
  }
  prop_suit
}

# Check if file exists
file_exist <- lapply(clim_scen_list, function(x) file.exists(x$path)) %>%
  unlist() %>%
  as.vector()

sum(file_exist) # How many files exist?

# Parallel processing to extract suitability scores at property centroids, 
# much quicker than polygon extraction for quick and dirty analysis
plan(multisession)
prop_suit_clim <- future.apply::future_lapply(clim_scen_list, 
                                              function(x) fcn_prop_suitability(x$path, gz = x$gz))
names(prop_suit_clim) <- lapply(clim_scen_list, function(x) x$id)
prop_suit_clim_df <- prop_suit_clim %>%
  bind_rows(.id = 'clim_scen') 

prop_suit_df <- cbind(NewPropID = properties[['NewPropID']], Area = properties[['ShapeArea']], prop_suit_clim_df)
write_csv(prop_suit_df, paste0(outpath, 'prop_suit_clim.csv'))
saveRDS(clim_scen_list, paste0(outpath, 'clim_scen_list.RData'))



