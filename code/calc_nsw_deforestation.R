# Calculate total woody vegetation extent in NSW

library(terra)
library(dplyr)
library(exactextractr)
library(sf)

wv <- terra::rast("data/s5hgps_nsw_y20082011_bcul0.tif")
nsw <- st_read("data/nsw_state_polygon_shp_gda2020/NSW_STATE_POLYGON_shp_GDA2020.shp") %>%
  st_transform(4283)
plot(nsw)


extract <- exact_extract(wv, nsw, 'sum')
