library(DescTools)

# Willingness-to-accept calculations
source('code/load_paths.R')

cost_pred <- read_csv("data/spatial_predictions_10yr.csv")

split_by_pct <- function(x, n_pct=5) {
  # Splits a random variable into a groups by percentile
  cutoffs <- quantile(x,0:(n_pct)/n_pct)
  return(cut(x, cutoffs, include.lowest = T, right = T))
}

set.seed(12305404)

rand_strat_sample <- function(cost_pred) {
  strata <- cost_pred %>%
    mutate(lval_perha_group = split_by_pct(LVAL/AREA, 20)) %>%
    group_by(KMR, lval_perha_group) %>%
    sample_frac(.05) %>%
    ungroup() %>%
    arrange(NewPropID)
  return(strata$NewPropID)
}

newpropid_rand <- replicate(100, rand_strat_sample(cost_pred)) %>%
  as.data.frame()
write_csv(newpropid_rand,"data/stratified_sample.csv")

nrow(strata)
sum(strata$AREA)

# 10-year costs (expressed in real units)
inflation_rate <- 0.02
inflation_vector <- (1+inflation_rate)^(seq(2025, 2085, 10) - 2020)
cost <- t(t(matrix(rep(properties_subset$MeanWTA * 1000 * (properties_subset$MeanProp * properties_subset$area_hectares), 7), ncol = 7)) * discount_vector) %>%
  as.data.frame()
colnames(cost) <- paste0("cost", as.character(seq(2025, 2085, 10)))
cost_df <- cbind(data.frame(NewPropID = properties_subset$NewPropID), cost)