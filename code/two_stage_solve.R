library(tidyverse)
library(stringr)
library(patchwork)
library(sf)
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
nsw_lga <- st_read( paste0(rdm_path, "\\planning_units\\planning_units.gdb"), layer = 'nsw_lga_pu_all') %>%
  st_transform(4326) %>%
  rename(lga = NSW_LGA__2)
properties_centroid <- properties %>%
  st_centroid() %>%
  st_transform(4326)
properties_intersect <- properties_centroid %>%
  st_join(nsw_lga)

selected_index <- clim_scen_list %>% lapply(function(x) x$id) %>% unlist() %>% as.vector()
# Select best 30 percent of sites by area based on average predictions
pct <- 0.3
properties_intersect <- cbind(prop_suit_df, 
                              lga = properties_intersect[['lga']], 
                              ha = properties_intersect[['Shape_Area.x']] / 10000)
properties_intersect$mean_suitability <- properties_intersect[selected_index] %>%
  rowMeans()
properties_intersect <- properties_intersect %>%
  group_by(lga) %>%
  arrange(-mean_suitability) %>%
  mutate(cum_ha = cumsum(ha)/sum(ha))

new_prop_id <- properties_intersect$NewPropID[!is.na(properties_intersect$lga) & properties_intersect$cum_ha < pct]
properties_intersect <- properties_intersect[properties_intersect$NewPropID %in% new_prop_id,]
prop_suit_df <- prop_suit_df[prop_suit_df$NewPropID %in% new_prop_id,] %>%
  arrange(NewPropID)
prop_suit_subset <- prop_suit_df[prop_suit_df$NewPropID %in% new_prop_id,selected_index]

# Weighted average of climate suitability
prop_ha <- properties[['Shape_Area']]/10000 # Area in hectares, to be replaced with koala habitat
prop_ha <- prop_ha[properties[['NewPropID']] %in% new_prop_id]

# 10-year costs
discount_factor <- 0.02
discount_vector <- (1+discount_factor)^(seq(2025, 2085, 10) - 2020)
prop_costs <- properties[properties$NewPropID %in% new_prop_id,] %>%
  arrange(NewPropID)
cost <- t(t(matrix(rep(prop_costs$MeanWTA, 7), ncol = 7)) * discount_vector) %>%
  as.data.frame()
colnames(cost) <- paste0("cost", as.character(seq(2025, 2085, 10)))
cost_df <- cbind(data.frame(NewPropID = prop_costs$NewPropID), cost)

# Suitability-weighted habitat
habitat_suitability_weighed <- prop_suit_df[selected_index] * prop_ha
ind <- 0
habitat_cc <- list()
new_prop_id_df <- data.frame(NewPropID = prop_costs$NewPropID)
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
write_csv(cost_df, "data/meanWTA_10yr.csv")
write_csv(habitat_cc_merged, "data/habitat_suitw_graham.csv")

# results <- fcn_two_stage_opt(cost = cost, quality = habitat_cc, area = prop_ha, t_prime = 4, target = rep(7000, 7),
#                              beta = 0, gamma = 0, type = 'recourse')





