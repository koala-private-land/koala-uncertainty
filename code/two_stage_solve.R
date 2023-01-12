library(tidyverse)
library(stringr)
library(patchwork)
library(sf)
library(mapview)

## Correlation between mean and variance
source('code/load_paths.R')
source('code/fcn_mv_opt.R')
prop_suit_df <- read_csv('data/prop_suit_clim.csv')
clim_scen_list <- readRDS('data/clim_scen_list.RData')

rcp <- lapply(clim_scen_list, function(x) x$rcp) %>% unlist() %>% as.vector()
year <- lapply(clim_scen_list, function(x) x$year) %>% unlist() %>% as.vector()
gcm <- lapply(clim_scen_list, function(x) x$gcm) %>% unlist() %>% as.vector()

# Load properties data
rdm_path <- "\\\\uq.edu.au\\UQ-Research\\KPRIVATE19-A2212\\data"
suitability_path <- "\\Graham_et_al_2019_koalas\\Phascolarctos_cinereus\\species_data_Phascolarctos_cinereus\\suitability\\"
properties_path <- "\\spatial_bid_model_predictions\\spatial_predictions_v1_0.gdb"
properties <- st_read(paste0(rdm_path, properties_path), layer = "spatial_pred_10yr_v1_0")
properties <- properties %>%
  mutate(area_meters = st_area(properties),
         area_hectares = units::set_units(area_meters, "hectares"),
         area_km = units::set_units(area_meters, "km^2"))
nsw_lga <- st_read( paste0(rdm_path, "\\planning_units\\planning_units.gdb"), layer = 'nsw_lga_pu_all') %>%
  st_transform(4326) %>%
  rename(lga = NSW_LGA__2)
properties_centroid <- properties %>%
  st_centroid() %>%
  st_transform(4326)
properties_intersect <- properties_centroid %>%
  st_join(nsw_lga)

selected_index <- clim_scen_list %>% lapply(function(x) x$id) %>% unlist() %>% as.vector()
# Select best 10 percent of sites by area based on average predictions
pct <- 0.1
properties_suit <- left_join(properties_intersect, prop_suit_df, by = 'NewPropID')
properties_suit <- properties_suit %>%
  mutate(mean_suitability = st_drop_geometry(properties_suit[selected_index]) %>% rowMeans())
properties_suit <- properties_suit %>%
  filter(MeanAdopt > 0.5) %>%
  group_by(lga) %>%
  arrange(-mean_suitability) %>%
  mutate(area_km = units::drop_units(area_km),
         cum_area_km = cumsum(area_km)/sum(area_km))

new_prop_id <- properties_suit$NewPropID[!is.na(properties_suit$lga) & properties_suit$cum_area_km < pct]
properties_subset <- properties_suit[properties_suit$NewPropID %in% new_prop_id,] %>%
  arrange(NewPropID)

# Extract koala habitat quantity 
khab_thres <- read_csv('data/khab_prop.csv') %>%
  rename(NewPropID = NEWPROPID) %>%
  mutate(khab_prop = SUM/COUNT)
properties_subset <- properties_subset %>%
  left_join(khab_thres, by = 'NewPropID') %>%
  mutate(khab_area_km = khab_prop * area_km)

# Weighted average of climate suitability
# 10-year costs (expressed in real units)
discount_factor <- 0.02
discount_vector <- (1+discount_factor)^(seq(2025, 2085, 10) - 2020)
cost <- t(t(matrix(rep(properties_subset$MeanWTA * 1000 * (properties_subset$MeanProp * properties_subset$area_hectares), 7), ncol = 7)) * discount_vector) %>%
  as.data.frame()
colnames(cost) <- paste0("cost", as.character(seq(2025, 2085, 10)))
cost_df <- cbind(data.frame(NewPropID = properties_subset$NewPropID), cost)

# Suitability-weighted habitat
habitat_suitability_weighed <- properties_subset[selected_index] * properties_subset$khab_area_km
ind <- 0
habitat_cc <- list()
new_prop_id_df <- data.frame(NewPropID = properties_subset$NewPropID)
for (g in unique(gcm)) {
  for (r in unique(rcp)) {
    ind <- ind+1
    model_index <- paste0(r, '_', g, '_', seq(2025,2085,10))
    habitat_cc[[ind]] <- cbind(new_prop_id_df, habitat_suitability_weighed[model_index])
    names(habitat_cc)[ind] <- paste0(r, '_', g)
    colnames(habitat_cc[[ind]]) <- c("NewPropID", seq(2025,2085,10))
  }
}
habitat_cc_merged <- bind_rows(habitat_cc, .id = 'climate_model')

## Export data to Julia to be solved by JuMP
saveRDS(properties_subset, file="data/properties_subset.RData")
write_csv(cost_df, "data/meanWTA_10yr.csv")
write_csv(habitat_cc_merged, "data/habitat_suitw_graham.csv")

# results <- fcn_two_stage_opt(cost = cost, quality = habitat_cc, area = prop_ha, t_prime = 4, target = rep(7000, 7),
#                              beta = 0, gamma = 0, type = 'recourse')




