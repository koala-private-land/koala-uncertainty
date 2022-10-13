library(nngeo)

## Make generic plots of the Graham et al SDMs
source('code/load_paths.R')
nsw_lga <- st_read( paste0(rdm_path, "\\planning_units\\planning_units.gdb"), layer = 'nsw_lga_pu_all') %>%
  st_union() %>%
  nngeo::st_remove_holes() %>%
  st_transform(4326) %>%
  st_union() %>%
  st_cast('POLYGON')

selected_scenarios <- c('RCP45_cccma-cgcm31_2015','RCP45_cccma-cgcm31_2085',
                        'RCP85_cccma-cgcm31_2015','RCP85_cccma-cgcm31_2085')
sdm <- lapply(clim_scen_list[selected_scenarios], function(x) {
  raster(x$path)
})

nsw_lga_sp <- sf::st_as_sf(nsw_lga)
sdm_masked <- lapply(sdm, function(x) raster::mask(x, nsw_lga_sp))

sdm_spdf <- lapply(sdm_masked, function(x) {
  df <- as(x, 'SpatialPixelsDataFrame') %>%
    as.data.frame()
  colnames(df) <- c('suitability', 'x', 'y')
  df
})
names(sdm_spdf) <- selected_scenarios
sdm_spdf_comb <- sdm_spdf %>%
  bind_rows(.id = 'clim_scen') %>%
  separate(clim_scen, c('rcp', 'model', 'year'), sep = '_') 

sdm_comb <- sdm_spdf_comb %>%
  ggplot(aes(x = x, y = y, fill = suitability)) +
  geom_raster() +
  facet_grid(rows = vars(rcp), cols = vars(year))+
  scale_fill_viridis_c(direction = -1) +
  coord_sf() +
  theme_void()
ggsave(sdm_comb, filename = 'plots/sdm_comb.png')  
