library(tidyverse)
library(stringr)
library(patchwork)
library(sf)
## Correlation between mean and variance
source('code/load_paths.R')
source('code/fcn_mv_opt.R')
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
properties_intersect <- properties_centroid %>%
  st_join(nsw_lga)

selected_index <- clim_scen_list[year == '2085'] %>% lapply(function(x) x$id) %>% unlist() %>% as.vector()
# Select best 10 percent of sites by area based on average predictions
pct <- 0.1
properties_intersect <- cbind(prop_suit_df, lga = properties_intersect[['lga']], ha = properties_intersect[['Shape_Area.x']] / 10000)
properties_intersect$mean_suitability <- properties_intersect[selected_index] %>%
  rowMeans()
properties_intersect <- properties_intersect %>%
  group_by(lga) %>%
  arrange(-mean_suitability) %>%
  mutate(cum_ha = cumsum(ha)/sum(ha))

new_prop_id <- properties_intersect$NewPropID[!is.na(properties_intersect$lga) & properties_intersect$cum_ha < pct]
properties_intersect <- properties_intersect[properties_intersect$NewPropID %in% new_prop_id,]
prop_suit_df <- prop_suit_df[prop_suit_df$NewPropID %in% new_prop_id,]
prop_suit_subset <- prop_suit_df[prop_suit_df$NewPropID %in% new_prop_id,selected_index]

# Weighted average of climate suitability
prop_ha <- properties[['Shape_Area']]/10000 # Area in hectares
prop_ha <- prop_ha[properties[['NewPropID']] %in% new_prop_id]

properties_intersect_lga <- properties_intersect %>%
  group_by(lga) %>%
  summarise_at(selected_index, funs(weighted.mean(., ha))) %>%
  arrange(lga) %>%
  dplyr::select('lga', selected_index)
nsw_prop_ha <- cbind(prop_suit_subset, lga = properties_intersect[['lga']], ha = prop_ha) %>%
  filter(!is.na(lga)) %>%
  group_by(lga) %>%
  summarise(ha = sum(ha))
lga_suit <- properties_intersect_lga
  
mean_sd_corr_df <- data.frame(lga = lga_suit[['lga']],
                              mean = rowMeans(lga_suit[selected_index]), 
                              sd = apply(lga_suit[selected_index], MARGIN=1, sd)) %>%
  left_join(nsw_prop_ha, by ='lga')
mean_sd_corr_lga <- nsw_lga %>%
  left_join(mean_sd_corr_df, by = 'lga') %>%
  st_simplify(dTolerance = 2000) %>%
  filter(!is.na(mean))

cutoff <- 8
mean_sd_corr_plot <- mean_sd_corr_lga %>%
  ggplot(aes(x = mean, y = sd)) +
  geom_point(aes(size = ha)) +
  ggrepel::geom_label_repel(aes(label = ifelse(mean/sd > cutoff, word(str_to_title(lga),1), NA))) +
  theme_bw()+
  scale_x_continuous('Habitat Suitability (Mean)') +
  scale_y_continuous('Habitat Suitability (SD)')
mean_sd_corr_plot

ggsave(mean_sd_corr_plot, filename = 'plots/mean_sd_corr_plot.png')


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
  ggrepel::geom_label_repel(aes(x = X, y = Y, label = ifelse(mean/sd > cutoff, word(str_to_title(lga), 1), NA)),nudge_x=-30) +
  scale_fill_viridis_c('Mean/SD', option = 'mako') +
  theme_void()
ggsave(mean_plot | sd_plot | sr_plot, filename = 'plots/lga_clim_risk.png', scale = 1.2)

## Mean Variance optimisation -------
lga_suit_ha <- lga_suit
for (s in selected_index) {
  lga_suit_ha[,s] <- lga_suit[,s] *mean_sd_corr_lga$ha / 100
}

lga_cov <- lga_suit_ha[selected_index] %>% t() %>% as.matrix() %>% cov()
target <- 0.1
obj <- -mean_sd_corr_lga$mean * mean_sd_corr_lga$ha
covmat <- as.matrix(lga_cov)
#A <- matrix(mean_sd_corr_lga$ha, nrow = 1, byrow = T)
#b <- c(sum(mean_sd_corr_lga$ha)*target)
A <- matrix(c(mean_sd_corr_lga$ha, -mean_sd_corr_lga$ha), nrow = 2, byrow = T)
b <- c(sum(mean_sd_corr_lga$ha)*target,-sum(mean_sd_corr_lga$ha)*target)
lambda_vec <- seq(0,1,0.1)
ev_results <- fcn_lp_opt(obj, A, b)
mv_results <- lapply(lambda_vec, function(l) fcn_mv_opt(obj, A, b, covmat, l))
names(mv_results) <- paste0('l', lambda_vec)
mv_decision <- lapply(mv_results, function(x) x$x) %>%
  bind_rows(.id = 'lambda')

mv_plot <- mean_sd_corr_lga %>%
  cbind(mv_decision) %>%
  pivot_longer(paste0('l', lambda_vec), names_to = 'lambda', values_to = 'decision') %>%
  mutate(lambda = factor(lambda, levels = paste0('l', lambda_vec), labels = paste0('λ = ', lambda_vec))) %>%
  ggplot(aes(fill = decision)) +
  geom_sf() +
  facet_wrap(~lambda,nrow = 2) +
  theme_void()
ggsave(mv_plot, filename = 'plots/lga_mv.png', scale = 1.5, units = 'mm', height = 120, width = 150)

benefit_dist <- t(as.matrix(lga_suit_ha[selected_index])) %*% as.matrix(mv_decision)
efficiency_frontier <- data.frame(mean = colMeans(benefit_dist), sd = apply(benefit_dist, 2, sd))
efficiency_frontier_plot <- ggplot(efficiency_frontier, aes(x = sd, y = mean)) +
  geom_point() +
  geom_smooth(se = F, color = '#4DBBD5FF') +
  scale_y_continuous('Expected benefits') +
  scale_x_continuous('Standard Deviation') +
  theme_bw()
ggsave(efficiency_frontier_plot, filename = 'plots/efficiency_frontier_plot.png', scale = 1, units = 'mm', height = 100, width = 120)
