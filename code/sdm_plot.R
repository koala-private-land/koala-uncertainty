library(nngeo)
library(raster)

## Make generic plots of the Graham et al SDMs
source('code/load_paths.R')
nsw_lga <- st_read( paste0(rdm_path, "\\planning_units\\planning_units.gdb"), layer = 'nsw_lga_pu_all') %>%
  st_union() %>%
  nngeo::st_remove_holes() %>%
  st_transform(4326) %>%
  st_union() %>%
  st_cast('POLYGON')

selected_scenarios <- c('RCP45_cccma-cgcm31_2015','RCP45_cccma-cgcm31_2045','RCP45_cccma-cgcm31_2085',
                        'RCP85_cccma-cgcm31_2015','RCP85_cccma-cgcm31_2045','RCP85_cccma-cgcm31_2085')
sdm <- lapply(clim_scen_list[selected_scenarios], function(x) {
  raster(x$path)
})

nsw_lga_sp <- sf::st_as_sf(nsw_lga)
sdm_masked <- lapply(sdm, function(x) raster::mask(x, nsw_lga_sp))

fcn_to_spdf <- function(e) {
  spdf <- lapply(e, function(x) {
  df <- as(x, 'SpatialPixelsDataFrame') %>%
    as.data.frame()
  colnames(df) <- c('suitability', 'x', 'y')
  df
  })
  names(spdf) <- selected_scenarios
  sdm_spdf_comb <- spdf %>%
    bind_rows(.id = 'clim_scen') %>%
    separate(clim_scen, c('rcp', 'model', 'year'), sep = '_') 
  }

sdm_masked_diff <- lapply(sdm_masked, function(x) x - sdm_masked[['RCP45_cccma-cgcm31_2015']])

sdm_spdf_comb <- fcn_to_spdf(sdm_masked)
sdm_spdf_diff <- fcn_to_spdf(sdm_masked_diff)

sdm_comb <- sdm_spdf_comb %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = suitability)) +
  geom_sf(data = nsw_lga, fill = NA) +
  facet_grid(rows = vars(rcp), cols = vars(year))+
  scale_fill_gradient(  low = "#ffffff", high = "#3C5488") +
  coord_sf() +
  theme_void()
ggsave(sdm_comb, filename = 'plots/sdm_comb.png')

sdm_comb_diff_plot <- ggplot() +
  geom_raster(data = sdm_spdf_diff, aes(x = x, y = y, fill = suitability)) +
  geom_sf(data = nsw_lga, fill = NA) +
  facet_grid(rows = vars(rcp), cols = vars(year))+
  scale_fill_gradient2('Change in \nsuitability index') +
  coord_sf() +
  theme_void()
sdm_comb_diff_plot
ggsave(sdm_comb_diff_plot, filename = 'plots/sdm_diff.png')
