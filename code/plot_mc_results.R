## Plots the results from sensitivity analyses

library(tidyverse)
library(ggplot2)
library(purrr)

results_dir <- "results/mc_sim/"

year_vec <- seq(2000, 2070, by=10)

# Baseline vs robust ------
default = c(1,4,0.25,1,10)
base_string <- "run_sp-%s_tt-%s_kt-%s_ns-%s_r-%s"
run_string <- do.call(sprintf, c(fmt = base_string, as.list(default)))
baseline_metric <- read_csv(paste0(results_dir, "metric_baseline_", run_string, ".csv"), col_types = cols())
default = c(1,4,0.25,12,10)
base_string <- "run_sp-%s_tt-%s_kt-%s_ns-%s_r-%s"
run_string <- do.call(sprintf, c(fmt = base_string, as.list(default)))
robust_metric <- read_csv(paste0(results_dir, "metric_nr_", run_string, ".csv"), col_types = cols())
flexible_metric <- read_csv(paste0(results_dir, 'metric_ar_', run_string, ".csv"), col_types = cols())

summary_stat <- function(metric, interval = c(0.05, 0.95)) {
  median <- apply(metric, 2, median)
  lb <- apply(metric, 2, quantile, interval[1])
  ub <- apply(metric, 2, quantile, interval[2])
  df <- data.frame(median, lb, ub, t = 1:ncol(metric))
  rownames(df) <-NULL
  return(df)
}
baseline_robust_df <- list(baseline = baseline_metric, robust = robust_metric, flexible = flexible_metric) %>%
  lapply(summary_stat) %>%
  bind_rows(.id = 'model') %>%
  mutate(year = year_vec[t]) %>%
  mutate(model = factor(model, c('baseline', 'robust', 'flexible'), c('SCP-CC', 'SCP-R', 'SCP-F'))) 

baseline_robust_df %>%
  ggplot(aes(x=year)) +
  geom_hline(yintercept = 7000) +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = model), alpha = 0.3) +
  geom_line(aes(y = median, color = model)) +
  geom_point(data = filter(baseline_robust_df, t == 3 & model == 'SCP-CC'), aes(y = median, color = model))+
  scale_y_continuous("Koala habitat (ha)") +
  scale_x_continuous("Year") +
  ggpubr::theme_pubr()
  


# Sensitivity plots ------
sp = 1 #1:10
tt = 4 # [1,2,3,4,5,6,7]
kt = 0.25 # [0.1, 0.15, 0.2, 0.25, 0.3]
ns = 1 # 1:12
rr = 50

plot_shaded_diff <- function(recourse = NULL, loop_vec = 1:10, default = c(1,4,0.25,1,10), dep_var = 1, interval = c(0.05, 0.95), realisation = NULL) {
  base_string <- "run_sp-%s_tt-%s_kt-%s_ns-%s_r-%s"
  cost_diff_pct <- lapply(loop_vec, function(x) {
    vec_var <- default
    vec_var[dep_var] <- x
    run_string <- do.call(sprintf, c(fmt = base_string, as.list(vec_var)))
    nr <- read_csv(paste0(results_dir, "cost_nr_", run_string, ".csv"), col_types = cols())
    ar <- read_csv(paste0(results_dir, "cost_ar_", run_string, ".csv"), col_types = cols()) # recourse solution
    tr <- read_csv(paste0(results_dir, "cost_tr_", run_string, ".csv"), col_types = cols()) # recourse solution
    fr <- read_csv(paste0(results_dir, "cost_fr_", run_string, ".csv"), col_types = cols()) # recourse solution
    diff <- list(ar = (ar - nr) / nr, tr = (tr - nr)/ nr, fr = (fr - nr)/nr)
  })
  names(cost_diff_pct) <- as.character(loop_vec)
  summary_diff_pct <- map(cost_diff_pct, function(list_x) {
    map(list_x, function(x) {
      col = x
      if (!is.null(realisation)) {
        col = x[,realisation]
      }
      df <- data.frame(median = median(unlist(col)), lb = quantile(unlist(col), interval[1]), ub = quantile(unlist(col), interval[2]))
      rownames(df) <- NULL
      df
    }) %>%
      bind_rows(.id = "recourse")
  }) %>%
    bind_rows(.id = 'name')
  
  summary_diff_pct %>%
    mutate(name = as.numeric(name)) %>%
    ggplot(aes(x = name, y = median)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = recourse), alpha = 0.25) +
    geom_hline(yintercept = 0) +
    geom_line(aes(color = recourse)) +
    scale_y_continuous("Change in conservation cost", labels = scales::unit_format(scale = 100,unit = "%")) +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    coord_cartesian(expand = F) +
    ggpubr::theme_pubr()
}

plot_shaded_diff(loop_vec = 1:10, dep_var = 1)
plot_shaded_diff(loop_vec = c(2,3,4,5,6,7), dep_var = 2, realisation = 1:10)
plot_shaded_diff(loop_vec = c(0.1, 0.15, 0.2, 0.25), dep_var = 3, realisation = 1:10)
plot_shaded_diff(loop_vec = 1:12, dep_var = 4, realisation = 1:10, interval = c(0.05, 0.95)) +
  scale_x_continuous("Number of climate scenarios")
