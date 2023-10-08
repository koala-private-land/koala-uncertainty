# Load shapes and save as shapes
rm(list=ls())
library(sf)

properties <- st_read("data/spatial_predictions_v1_1.gdb", layer = "spatial_pred_inperp_v1_1")
prop_centroid <- properties %>%
  st_centroid() %>%
  select(NewPropID, MeanProp, UpperProp, Shape_Area) %>%
  mutate(AREA = Shape_Area * 0.0001)
prop_df <- st_drop_geometry(prop_centroid) %>%
  cbind(st_coordinates(prop_centroid))
## Spatial plots
aus_border <- geodata::gadm(country = "AUS", level = 1, path = "data/", resolution = 2) %>%
  sf::st_as_sf()
aus_border_union <- aus_border %>% st_set_precision(1e6) %>% st_make_valid() %>% st_union() %>% plot()
nsw_lga <- st_read( "data/planning_units.gdb", layer = 'nsw_lga_pu_all')
nsw_lga_union <- nsw_lga%>%
  st_union() %>%
  nngeo::st_remove_holes() %>%
  sf::st_transform(4326)
nsw_bbox <- st_bbox(nsw_lga_union)
bbox_buffer <- 1
aus_xlim <- c(113.338953078, 153.569469029)
aus_ylim <- c(-43.6345972634,-10.6681857235)
prop_lga <- st_join(prop_centroid, nsw_lga)
prop_lga_lookup <- prop_lga %>% 
  st_drop_geometry() %>%
  select(NewPropID, NSW_LGA__2)

save.image('plots/shapes.RData')
rm(list=ls())
gc()