library(stringr)
library(patchwork)

## Correlation between mean and variance
prop_suit_df <- read_csv('../../Koala/temporal_planning/prop_suit_clim.csv')
clim_scen_list <- readRDS('../../Koala/temporal_planning/clim_scen_list.RData')

rcp <- lapply(clim_scen_list, function(x) x$rcp) %>% unlist() %>% as.vector()
year <- lapply(clim_scen_list, function(x) x$year) %>% unlist() %>% as.vector()
gcm <- lapply(clim_scen_list, function(x) x$gcm) %>% unlist() %>% as.vector()

# Load properties data
rdm_path <- "\\\\uq.edu.au\\UQ-Research\\KPRIVATE19-A2212\\data"
suitability_path <- "\\Graham_et_al_2019_koalas\\Phascolarctos_cinereus\\species_data_Phascolarctos_cinereus\\suitability\\"
properties_path <- "\\spatial_bid_model_predictions\\spatial_predictions_v1_0.gdb"
properties <- st_read(paste0(rdm_path, properties_path), layer = "spatial_pred_inperp_v1_0")
nsw_lga <- st_read( paste0(rdm_path, "\\planning_units\\planning_units.gdb"), layer = 'nsw_lga_pu_all') %>%
  st_transform(4326) %>%
  rename(lga = NSW_LGA__2)
properties_centroid <- properties %>%
  st_centroid() %>%
  st_transform(4326)
properties_intersect <- properties_centroid  %>%
  st_join(nsw_lga)

selected_index <- clim_scen_list[year == '2085' & rcp == 'RCP45'] %>% lapply(function(x) x$id) %>% unlist() %>% as.vector()
prop_suit_subset <- prop_suit_df[selected_index]

# Weighted average of climate suitability
prop_ha <- properties[['Shape_Area']]/10000 # Area in hectares
lga_suit <- cbind(prop_suit_subset, lga = properties_intersect[['lga']], ha = prop_ha) %>%
  group_by(lga) %>%
  summarise_at(selected_index, funs(weighted.mean(., ha))) %>%
  arrange(lga) %>%
  dplyr::select('lga', selected_index)
mean_sd_corr_df <- data.frame(lga = lga_suit[['lga']],
                              mean = rowMeans(lga_suit[selected_index]), 
                              sd = apply(lga_suit[selected_index], MARGIN=1, sd))
mean_sd_corr_lga <- nsw_lga %>%
  left_join(mean_sd_corr_df, by = 'lga') %>%
  mutate(ha = Shape_Area / 10000) %>%
  st_simplify(dTolerance = 2000)

mean_sd_corr_lga %>%
  ggplot(aes(x = mean, y = sd)) +
  geom_point(aes(size = ha)) +
  ggrepel::geom_label_repel(aes(label = ifelse(mean/sd > 12, word(str_to_title(lga), 1), NA))) +
  theme_bw()+
  scale_x_continuous('Habitat Suitability (Mean)') +
  scale_y_continuous('Habitat Suitability (SD)')

mean_sd_corr_lga[c('X', 'Y')] <- mean_sd_corr_lga %>% st_centroid() %>% st_coordinates()

mean_plot <- mean_sd_corr_lga %>%
  ggplot() +
  geom_sf(aes(fill = mean)) +
  #ggrepel::geom_label_repel(aes(x = X, y = Y, label = ifelse(mean/sd > 12, word(str_to_title(lga), 1), NA)),nudge_x=-30) +
  scale_fill_viridis_c('Mean', direction = -1, option = 'viridis') +
  theme_void()
sd_plot <- mean_sd_corr_lga %>%
  ggplot() +
  geom_sf(aes(fill = sd)) +
  #ggrepel::geom_label_repel(aes(x = X, y = Y, label = ifelse(mean/sd > 12, word(str_to_title(lga), 1), NA)),nudge_x=-30) +
  scale_fill_viridis_c('SD', direction = -1, option = 'magma') +
  theme_void()
sr_plot <- mean_sd_corr_lga %>%
  ggplot() +
  geom_sf(aes(fill = mean/sd)) +
  ggrepel::geom_label_repel(aes(x = X, y = Y, color = 'white', label = ifelse(mean/sd > 12, word(str_to_title(lga), 1), NA)),nudge_x=-30) +
  scale_fill_viridis_c('Mean/SD', direction = -1, option = 'mako') +
  theme_void()
ggsave(mean_plot | sd_plot | sr_plot, filename = 'plots/lga_clim_risk.png', scale = 1.2)


mean_suit <- rowMeans(prop_suit_subset)
sd_suit <- apply(prop_suit_subset, MARGIN = 1, sd)
mean_sd_corr_df <- cbind(mean = mean_suit, sd = sd_suit, prop_size = prop_size,
                         nsw_lga = properties_intersect[['NSW_LGA__2']]) %>%
  as.data.frame()
mean_sd_corr_lga <- mean_sd_corr_df %>%
  summarise()
  
mean_sd_corr_df %>%
  filter(prop_size > 174961) %>%
  ggplot(aes(x = mean, y = sd)) +
  geom_point() +
  theme_bw()

