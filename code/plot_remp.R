library(terra)
library(tidyterra)
library(ggplot2)
library(sf)

load('plots/shapes.RData')

remp_model_path <- "M:\\Users\\uqfcho\\Documents\\remp_models_2023\\data\\Combine_Inland_Coastal"
models <- c('CCCMA_R1', "CCCMA_R2", "CCCMA_R3", "CSIRO_R1", "CSIRO_R2", "CSIRO_R3",
            "ECHAM_R1", "ECHAM_R2", "ECHAM_R3", "MIROC_R1", "MIROC_R2", "MIROC_R3")
#models <- models[1]
tif_paths_2070 <- sapply(models, function(m) {
  pattern <- paste0(remp_model_path, '\\%s\\Occupancy\\Koala_%s_Stable_Coupled_Combined_Pi_t6.tif')
  path <- sprintf(pattern, m,m)
})
tif_2070 <-  terra::rast(tif_paths_2070)
names(tif_2070) <- models
tif_agg_2070 <- terra::aggregate(tif_2070, fact = 20, fun = "mean", na.rm = TRUE)

threshold <- 0.25

classify_matrix <- c(0, threshold, 0,
                     threshold, 1, 1) %>%
  matrix(byrow = T, ncol = 3)

nsw_lga_trans <- terra::vect(nsw_lga_union) %>% terra::project(tif_agg_2070)

remp_threshold <- terra::classify(tif_agg_2070, classify_matrix)
remp_threshold <- terra::mask(remp_threshold,nsw_lga_trans)
remp_threshold <- terra::crop(remp_threshold,nsw_lga_trans)

remp_plot <- ggplot() +
  geom_spatraster(data=remp_threshold, na.rm = TRUE) +
  geom_spatvector(data = nsw_lga_trans, na.rm = TRUE, linewidth = 0.5, fill = NA) +
  scale_fill_whitebox_b(breaks = c(0, 0.5, 1), palette = "atlas", direction = -1)+
  facet_wrap(~lyr) +
  theme_void()

ggsave("plots/remp_plot.png", remp_plot, width = 3000, height = 3000, units = 'px')
