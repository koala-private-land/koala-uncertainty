library(dplyr)
library(readr)

# Willingness-to-accept calculations
source('code/load_paths.R')

cost_pred <- read_csv("data/spatial_predictions_inf.csv")

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
    sample_frac(.25) %>%
    ungroup() %>%
    arrange(NewPropID)
  return(strata$NewPropID)
}

newpropid_rand <- replicate(100, rand_strat_sample(cost_pred)) %>%
  as.data.frame()
write_csv(newpropid_rand,"data/stratified_sample.csv")