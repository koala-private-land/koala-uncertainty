## Plots the results from sensitivity analyses

library(tidyverse)
library(ggplot2)
library(purrr)
library(patchwork)
library(geodata)
library(sf)
library(ggthemes)
library(ggplotify)
library(eulerr)

results_dir <- "results/mc_sim_mcmc/"
#results_dir <- "~/Documents/OneDrive - The University of Queensland/Documents/GitHub/koala-uncertainty/results/mc_sim/"

# Default parameters
sp = 1 #1:10
tt = 6 # [0,1,2,3,4,5,6,7]
kt = 0.25 # [0.1, 0.15, 0.2, 0.25, 0.3]
ns = 12 # 1:12
rr = 30
sdr = 0.02
dr = 0.1
k = 7000

colorpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442" ,"#0072B2")

# Fill
scale_fill_colorblind7 = function(...){
  scale_fill_discrete(..., type = colorpal)
}

# Color
scale_color_colorblind7 = function(...){
  scale_color_discrete(..., type = colorpal)
}

default_params = c(k, sp, tt, kt, ns, rr, sdr, dr)
base_string <- "run_k-%s_sp-%s_tt-%s_kt-%s_ns-%s_r-%s_sdr-%s_dr-%s"

flexibility_differences <- function(recourse = c(T,T,F,F), loop_vec = 1:10, params = default_params, 
                                    dep_var = 1, interval = c(0.05, 0.95), realisation = NULL, base_dir = results_dir) {
  # dep_var starts at 0
  cost_diff_pct <- lapply(loop_vec, function(x) {
    vec_var <- params
    vec_var[dep_var] <- x
    run_string <- do.call(sprintf, c(fmt = base_string, as.list(vec_var)))
    diff <- list()
    if (recourse[1]) {
      nr <- read_csv(paste0(results_dir, "cost_nr_", run_string, ".csv"), col_types = cols())
    }
    if (recourse[2]){
      ar <- read_csv(paste0(results_dir, "cost_ar_", run_string, ".csv"), col_types = cols()) # recourse solution
      diff$ar = (ar - nr) / nr
    }
    if (recourse[3]){
      tr <- read_csv(paste0(results_dir, "cost_tr_", run_string, ".csv"), col_types = cols()) # recourse solution
      diff$tr = (tr - nr)/ nr
      }
    if (recourse[4]) {
      fr <- read_csv(paste0(results_dir, "cost_fr_", run_string, ".csv"), col_types = cols()) # recourse solution
      diff$fr = (fr - nr)/nr
    }
    diff
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
  
  return(summary_diff_pct)
}

year_vec <- seq(2000, 2070, by=10)

# Baseline vs robust ------

# Load property data
load('plots/shapes.RData')

aus_plot <- ggplot() +
  geom_sf(data = aus_border, fill = 'gray80', color = 'white') +
  geom_sf(data = nsw_lga_union, size = 0.5, color = NA, fill = '#FFD580') +
  geom_rect(aes(xmin = nsw_bbox$xmin - bbox_buffer, 
                xmax = nsw_bbox$xmax + bbox_buffer, 
                ymin = nsw_bbox$ymin - bbox_buffer, 
                ymax = nsw_bbox$ymax + bbox_buffer), 
            fill = NA, color = 'gray50') +
  coord_sf(xlim = aus_xlim, ylim = aus_ylim) +
  theme_void()
#theme_bw() +
#theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
#      panel.background = element_rect(fill = 'aliceblue'))

spatial_pred <- read_csv('data/spatial_predictions_10yr.csv')

scen_list <- c('Inflexible - Ignore Risk', 'Inflexible - Robust', 'Flexible - Climate Change', 'Flexible & Learning - Climate Change')
plot_list <- list()

fcn_decision_set <- function(a,b,area) {
  # Set of decisions
  x_str <- paste0("x", 1:rr)
  ax <- a[x_str]
  bx <- b[x_str]
  u_idx <- ax * bx # union
  u <- colSums(u_idx * area)
  ao <- colSums(ax * (u_idx<0.01) * area)
  bo <- colSums(bx * (u_idx<0.01) * area)
  area_sum <- u+ao+bo
  return(data.frame(u=u,a=ao,b=bo))
}

scen_color_def <- as.vector(colorpal[1:length(scen_list)])
names(scen_color_def) <- scen_list

for (i in 4:length(scen_list)) {
  scen_list_i <- scen_list[1:i]

  get_run_string <- function(param=default) do.call(sprintf, c(fmt = base_string, as.list(param)))
  baseline_string <- get_run_string(c(k, sp, tt, kt, ns, rr, sdr, dr))
  robust_string   <- get_run_string(c(k, sp, tt, kt, ns, rr, sdr, dr))
  flexible_string <- get_run_string(c(k, sp, tt, kt, ns, rr, sdr, dr))
  flexible_learning_string <- get_run_string(c(k, sp, tt, kt, 1, rr, sdr, dr))
  
  baseline_metric <- read_csv(paste0(results_dir, "metric_baseline_", baseline_string, ".csv"), col_types = cols())
  robust_metric <- read_csv(paste0(results_dir, "metric_nr_", robust_string, ".csv"), col_types = cols())
  flexible_metric <- read_csv(paste0(results_dir, 'metric_ar_', flexible_string, ".csv"), col_types = cols())
  flexible_learning_metric <- read_csv(paste0(results_dir, 'metric_ar_', flexible_learning_string, ".csv"), col_types = cols())
  
  summary_stat <- function(metric, interval = c(0.05, 0.95)) {
    median <- apply(metric, 2, median)
    lb <- apply(metric, 2, quantile, interval[1])
    ub <- apply(metric, 2, quantile, interval[2])
    df <- data.frame(median, lb, ub, t = 1:ncol(metric))
    rownames(df) <-NULL
    return(df)
  }
  
  baseline_robust_df <- list(baseline = baseline_metric, robust = robust_metric, flexible = flexible_metric, flexible_learning = flexible_learning_metric) %>%
    lapply(summary_stat) %>%
    bind_rows(.id = 'model') %>%
    mutate(year = year_vec[t]) %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), scen_list)) %>%
    filter(model %in% scen_list_i)
  
  year_trend_plot <- baseline_robust_df %>%
    filter(year >= 2020) %>%
    filter(model %in% scen_list_i) %>%
    ggplot(aes(x=year)) +
    geom_hline(yintercept = 7000) +
    geom_vline(xintercept = 2020, linetype = 2, color = 'gray50') +
    geom_vline(xintercept = year_vec[tt]-5, linetype = 2, color = 'gray50') +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = model), alpha = 0.25) +
    geom_line(aes(y = median, color = model), linewidth = 1) +
    geom_point(data = filter(baseline_robust_df, t %in% c(3)), aes(y = median, color = model), size = 2)+
    geom_point(data = filter(baseline_robust_df, t %in% tt & model %in% scen_list[3:4]), aes(y = median, color = model), size = 2)+
    scale_y_continuous("Koala habitat (ha)") +
    scale_x_continuous("Year") +
    scale_color_colorblind7() +
    scale_fill_colorblind7() +
    coord_cartesian(xlim = c(2020,2070), ylim = c(0, 18000)) +
    guides(color = 'none', fill = 'none')+
    ggpubr::theme_pubr()
  
  year_trend_split_plot <- baseline_robust_df %>%
    filter(year >= 2020) %>%
    filter(model %in% scen_list_i) %>%
    ggplot(aes(x=year)) +
    geom_hline(yintercept = 7000) +
    geom_vline(xintercept = 2020, linetype = 2, color = 'gray50') +
    geom_vline(xintercept = year_vec[tt]-5, linetype = 2, color = 'gray50') +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = model), alpha = 0.25) +
    geom_line(aes(y = median, color = model), linewidth = 1) +
    geom_point(data = filter(baseline_robust_df, t %in% c(3)), aes(y = median, color = model), size = 2)+
    geom_point(data = filter(baseline_robust_df, t %in% tt & model %in% scen_list[3:4]), aes(y = median, color = model), size = 2)+
    scale_y_continuous("Koala habitat (ha)") +
    scale_x_continuous("Year") +
    scale_color_colorblind7() +
    scale_fill_colorblind7() +
    facet_grid(~model, labeller = label_wrap_gen()) +
    coord_cartesian(xlim = c(2020,2070), ylim = c(0, 18000)) +
    guides(color = 'none', fill = 'none')+
    ggpubr::theme_pubr()
  
  end_range_plot <- baseline_robust_df %>%
    filter(year==2070) %>%
    #group_by(model) %>%
    #summarise(median = median(median), lb = min(lb), ub = max(ub)) %>%
    mutate(name_num = as.numeric(as.factor(model))) %>%
    filter(model %in% scen_list_i) %>%
    ggplot(aes(x = name_num, y = median, fill = model)) +
    geom_hline(yintercept = 7000) +
    geom_rect(aes(ymin = lb, ymax = ub, xmin = name_num - 0.4, xmax = name_num + 0.4), color = 'gray80', linewidth = 0.5)+
    geom_segment(aes(y = median, yend = median, x = name_num-0.4, xend = name_num+0.4), color = 'white', linewidth = 0.5) +
    #geom_text(aes(y = 0, x = name_num, label = model, color = model), size = 4, vjust = 0.5) +
    coord_cartesian(ylim = c(0, 18000)) +
    guides(color = 'none') +
    scale_fill_manual("Strategy",values = scen_color_def, labels = function(x) str_wrap(x, width = 24)) +
    scale_color_manual("Strategy",values = scen_color_def, labels = function(x) str_wrap(x, width = 24)) +
    theme_void() +
    theme(axis.line.x = element_blank())
  
  ## Import decision vectors
  decision_vector_calc <- function(df) {
    cutoff <- rr * 0.1 # Proportion chosen across simulations to be in the plot
    sum_func <- function(v) mean(v)
    df$x <- apply(df[,paste0('x', 1:rr)], 1, sum_func)
    df$y <- apply(df[,paste0('sum_y', 1:rr)] / 12, 1, sum_func)
    df$w <- apply(df[,paste0('sum_w', 1:rr)] / 12, 1, sum_func)
    df
  }
  
  baseline_decisions <- read_csv(paste0(results_dir, "decision_baseline_", baseline_string, ".csv"), col_types = cols()) %>% 
    decision_vector_calc()
  robust_decisions <- read_csv(paste0(results_dir, "decision_nr_", robust_string, ".csv"), col_types = cols()) %>% 
    decision_vector_calc()
  flexible_decisions <- read_csv(paste0(results_dir, 'decision_ar_', flexible_string, ".csv"), col_types = cols()) %>% 
    decision_vector_calc()
  flexible_learning_decisions <- read_csv(paste0(results_dir, 'decision_ar_', flexible_learning_string, ".csv"), col_types = cols()) %>% 
    decision_vector_calc()
  
  ## Total conservation costs
  baseline_costs <- read_csv(paste0(results_dir, "cost_baseline_", baseline_string, ".csv"), col_types = cols())
  robust_costs <- read_csv(paste0(results_dir, "cost_nr_", baseline_string, ".csv"), col_types = cols())
  flexible_costs <- read_csv(paste0(results_dir, 'cost_ar_', flexible_string, ".csv"), col_types = cols())
  flexible_learning_costs <- read_csv(paste0(results_dir, 'cost_ar_', flexible_learning_string, ".csv"), col_types = cols())
  
  ## Area of properties
  baseline_area <- read_csv(paste0(results_dir, "area_", baseline_string, ".csv"), col_types = cols()) 
  robust_area <- read_csv(paste0(results_dir, "area_", baseline_string, ".csv"), col_types = cols())
  flexible_area <- read_csv(paste0(results_dir, "area_", baseline_string, ".csv"), col_types = cols()) 
  flexible_learning_area <- read_csv(paste0(results_dir, "area_", baseline_string, ".csv"), col_types = cols())
  
  ## Upper limit of area on offer (offered area)
  baseline_area_max <- apply(baseline_area, 1, max)
  robust_area_max <- apply(robust_area, 1, max)
  flexible_area_max <- apply(flexible_area, 1, max)
  flexible_learning_area_max <- apply(flexible_learning_area, 1, max)
  
  ## Sum of area
  calculate_offered_area <- function(decisions, area) {
    return(as.data.frame(t(apply(area, 1, max)) %*% as.matrix(decisions[,c('x','y','w')])))
  }
  calculate_area <- function(area, decisions) {
    x <- decisions[,paste0('x', 1:rr)]
    y <- decisions[,paste0('sum_y', 1:rr)]/12
    w <- decisions[,paste0('sum_w', 1:rr)]/12
    
    data.frame(x = (area * x) %>% colSums() %>% as.vector(),
               y = (area * y) %>% colSums() %>% as.vector(),
               w = (area * w) %>% colSums() %>% as.vector())
  }
  
  area_range <- function(area) {
    area <- as.data.frame(area)
    median = apply(area, 2, quantile, probs = .5)
    lb = apply(area, 2, quantile, probs = .05)
    ub = apply(area, 2, quantile, probs = .95)
    df <- rbind(median = median, lb = lb, ub = ub) %>% as.data.frame()
    df$stat <- as.list(row.names(df))
    return(df)
  }
  
  baseline_sum_area <- calculate_area(baseline_area, baseline_decisions) %>% area_range()
  robust_sum_area <- calculate_area(robust_area, robust_decisions) %>% area_range()
  flexible_sum_area <- calculate_area(flexible_area, flexible_decisions) %>% area_range()
  flexible_learning_sum_area <- calculate_area(flexible_learning_area, flexible_learning_decisions) %>% area_range()
  
  area_to_ts <- function(area) {
    area <- as.data.frame(area)
    year_vec <- seq(2000, 2070, by=10)
    total_area <- list()
    for (i in 1:nrow(area)) {
      total_area[[i]] <- data.frame(year = year_vec, area = rep(area$x[i], length(year_vec)) + c(rep(0, tt-1), rep(area$y[i] - area$w[i], length(year_vec) - tt + 1)))
    }
    names(total_area) <- area$stat
    return(bind_rows(total_area, .id = 'stat'))
  }
  
  baseline_area_ts <- area_to_ts(baseline_sum_area)
  robust_area_ts <- area_to_ts(robust_sum_area)
  flexible_area_ts <- area_to_ts(flexible_sum_area)
  flexible_learning_area_ts <- area_to_ts(flexible_learning_sum_area)
  
  baseline_offered_area <- calculate_offered_area(baseline_decisions, baseline_area) %>% area_to_ts()
  robust_offered_area <- calculate_offered_area(robust_decisions, robust_area) %>% area_to_ts()
  flexible_offered_area <- calculate_offered_area(flexible_decisions, flexible_area)%>% area_to_ts()
  flexible_learning_offered_area <- calculate_offered_area(flexible_learning_decisions, flexible_learning_area) %>% area_to_ts()
  
  offered_area <- list(baseline = baseline_offered_area, robust = robust_offered_area, 
                       flexible = flexible_offered_area, flexible_learning = flexible_learning_offered_area) %>%
    bind_rows(.id = 'model') %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), scen_list)) %>%
    filter(year >= 2020) %>%
    filter(model %in% scen_list_i)
  
  covenanted_area_plot <- list(baseline = baseline_area_ts, robust = robust_area_ts, flexible = flexible_area_ts, flexible_learning = flexible_learning_area_ts) %>%
    bind_rows(.id = 'model') %>%
    pivot_wider(names_from = 'stat', values_from = 'area') %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), scen_list)) %>%
    filter(year >= 2020) %>%
    filter(model %in% scen_list_i) %>%
    ggplot(aes(x = year)) +
    geom_vline(xintercept = 2020, linetype = 2, color = 'gray50') +
    geom_vline(xintercept = year_vec[tt]-5, linetype = 2, color = 'gray50') +
    geom_ribbon(aes(fill = model, ymin = lb, ymax = ub), alpha = 0.25) +
    geom_line(aes(color = model, y = median), size = 1) +
    #geom_line(aes(color = model, y = median)) +
    scale_y_continuous("Covenant area (ha)") +
    scale_x_continuous("Year") +
    scale_color_colorblind7() +
    scale_fill_colorblind7() +
    coord_cartesian(xlim = c(2020,2070), ylim = c(3000, 15000)) +
    guides(color = 'none', fill = 'none')+
    ggpubr::theme_pubr() +
    annotate('text', x = mean(c(year_vec[3], year_vec[tt]-5)), y = 15000, color = 'black', label = 'Stage 1', vjust=1) +
    annotate('text', x = mean(c(year_vec[tt]-5, year_vec[length(year_vec)])) , y = 15000, color = 'black', label = 'Stage 2', vjust=1) +
    theme(axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  offered_area_end <- offered_area %>%
    mutate(name_num = as.numeric(model)) %>%
    filter(year == 2070)
  
  covenanted_area_stages_full <- list(baseline = baseline_area_ts, robust = robust_area_ts, flexible = flexible_area_ts, flexible_learning = flexible_learning_area_ts) %>%
    bind_rows(.id = 'model') %>%
    pivot_wider(names_from = 'stat', values_from = 'area') %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), scen_list)) 
  
  covenanted_area_stages <- covenanted_area_stages_full %>%
    filter(year %in% c(2020,2070)) %>%
    filter(model %in% scen_list_i) %>%
    mutate(name_num = as.numeric(as.factor(model))) %>%
    mutate(stage = factor(year, c(2020,2070), c('Stage 1', 'Stage 2'))) %>%
    mutate(stage_num = ifelse(stage == 'Stage 1', 1, 2)) %>%
    mutate(xmin = stage_num + (name_num-2.5)*.2- 0.05, xmax = stage_num + (name_num-2.5)*.2 + 0.05) 
  
  covenanted_area_ref <- covenanted_area_stages %>%
    filter(stage == 'Stage 1' & model %in% scen_list[3:4]) %>%
    mutate(xmin = xmin + 1, xmax = xmax + 1)
  
  covenanted_area_ref2 <- covenanted_area_stages %>%
    filter(model %in% scen_list[3:4]) %>%
    mutate(xmin = xmin + 1, xmax = xmax + 1)
  
  covenanted_area_stages_plot <- covenanted_area_stages %>%
    ggplot(aes(x = stage_num, y = median, color = model)) +
    geom_rect(aes(fill = model, ymin = lb, ymax = ub, xmin = xmin, xmax = xmax)) +
    geom_segment(aes(y = median, yend = median, x = xmin-0.01, xend = xmax+0.01), color = 'white', size = 1) +
    geom_segment(data = covenanted_area_ref, aes(y = median, yend = median, x = xmin-0.01, xend = xmax+0.075, color = model)) +
    annotate('segment', y = filter(covenanted_area_ref2, stage == "Stage 1")$median, yend = filter(covenanted_area_ref2, stage == "Stage 2")$median, 
             x = covenanted_area_ref$xmax + 0.04, xend = covenanted_area_ref$xmax + 0.04, 
             color = colorpal[3:4],
             size = .75,
             arrow = arrow(length = unit(0.2, "cm"))) +
    geom_vline(xintercept = c(0.5,1.5), linetype = 2, color = 'gray50') +
    scale_fill_colorblind7() +
    scale_color_colorblind7() +
    guides(color = 'none', fill = 'none') +
    scale_y_continuous("Covenant area (ha)") +
    ggpubr::theme_pubr() +
    annotate('text', x = 1, y = 18000, color = 'black', label = 'Stage 1', vjust=1) +
    annotate('text', x = 2 , y = 18000, color = 'black', label = 'Stage 2', vjust=1) +
    coord_cartesian(xlim = c(0.5,2.5), ylim = c(6000,18000)) +
    theme(axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) 
  
  covenanted_area_end_plot <- list(baseline = baseline_area_ts, robust = robust_area_ts, flexible = flexible_area_ts, flexible_learning = flexible_learning_area_ts) %>%
    bind_rows(.id = 'model') %>%
    pivot_wider(names_from = 'stat', values_from = 'area') %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), scen_list)) %>%
    filter(year == 2070) %>%
    filter(model %in% scen_list_i) %>%
    mutate(name_num = as.numeric(as.factor(model))) %>%
    ggplot(aes(x = name_num, y = median, fill = model)) +
    #geom_hline(yintercept = 7000) +
    geom_rect(aes(ymin = lb, ymax = ub, xmin = name_num - 0.4, xmax = name_num + 0.4))+
    geom_segment(aes(y = median, yend = median, x = name_num-0.4, xend = name_num+0.4), color = 'white', size = 1) +
    #geom_point(data = offered_area_end, aes(y = area, x = name_num, color = model), size = 2, shape = 23) +
    #geom_text(aes(y = 3500, x = name_num, label = model, color = model), size = 4, vjust = 0) +
    coord_cartesian(ylim = c(3000, 15000)) +
    guides(color = 'none', fill = 'none') +
    scale_fill_colorblind7() +
    scale_color_colorblind7() +
    theme_void() +
    theme(axis.line.x = element_blank())
  
  year_split_inset <- baseline_robust_df %>%
    filter(year >= 2020) %>%
    ggplot(aes(x=year)) +
    geom_hline(yintercept = 7000) +
    geom_vline(xintercept = 2020, linetype = 1, color = 'gray50') +
    geom_vline(xintercept = year_vec[tt]-5, linetype = 1, color = 'gray50') +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = model), alpha = 0.5) +
    geom_line(data = filter(covenanted_area_stages_full, year >=2020), 
              aes(x = year, y = lb, color = model), linetype = 3, linewidth = 0.5)+
    geom_line(data = filter(covenanted_area_stages_full, year >=2020), 
              aes(x = year, y = ub, color = model), linetype = 3, linewidth = 0.5)+
    geom_line(aes(y = median, color = model), linewidth = 1) +
    geom_point(data = filter(baseline_robust_df, t %in% c(3)), aes(y = median, color = model), size = 2)+
    geom_point(data = filter(baseline_robust_df, t %in% tt & model %in% scen_list[3:4]), aes(y = median, color = model), size = 2)+
    annotate("text", label = "Stage 1", x = mean(c(2020,year_vec[tt]-5)), y = 2000+16500, hjust = 0.5, vjust = 1, color = 'gray50') +
    annotate("text", label = "Stage 2", x = mean(c(2070,year_vec[tt]-5)), y = 2000+16500, hjust = 0.5, vjust = 1, color = 'gray50') +
    annotate("segment", x = 2020, xend = year_vec[tt]-5, y = 200+16500, yend = 200+16500, arrow = arrow(ends = "both",length = unit(.2,"cm")), color = 'gray50') +
    annotate("segment", x = 2070, xend = year_vec[tt]-5, y = 200+16500, yend = 200+16500, arrow = arrow(ends = "both",length = unit(.2,"cm")), color = 'gray50') +
    scale_y_continuous("Policy target", labels = function(x) paste0(x*100 / 7000, '%'), breaks = ((0:8)*50) * 70) +
    scale_x_continuous("Year", breaks = c(2020, 2040, 2060)) +
    scale_color_manual(values = scen_color_def) +
    scale_fill_manual(values = scen_color_def) +
    facet_grid(~model, labeller = label_wrap_gen()) +
    coord_cartesian(xlim = c(2020,2070), ylim = c(0, 18000)) +
    guides(color = 'none', fill = 'none') +
    ggpubr::theme_pubr() +
    theme(strip.background = element_blank())
  
  ## Cost plot ----------
  cost_list <- list(baseline = baseline_costs, robust = robust_costs, flexible = flexible_costs, flexible_learning = flexible_learning_costs)
  costs <- cost_list%>%
    lapply(function(x, interval = c(0.05, 0.95)) {
      x_unlist <- unlist(x)
      median <- median(x_unlist)
      lb <- quantile(x_unlist, interval[1])
      ub <- quantile(x_unlist, interval[2])
      df <- data.frame(stat = c('median', 'lb', 'ub'), cost = c(median, lb, ub))
      row.names(df) <- NULL
      return(df)
    }) %>%
    bind_rows(.id = 'model') %>%
    pivot_wider(names_from = 'stat', values_from = 'cost')
  
  cost_dist <- cost_list %>%
    lapply(function(x) unlist(x)) %>%
    bind_cols(.id = 'model') %>%
    pivot_longer(all_of(names(cost_list)), names_to = 'model', values_to = 'cost') %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), scen_list)) %>%
    ggplot(aes(fill = model, y = model, x = cost)) +
    ggridges::geom_density_ridges(color = 'white') +
    scale_fill_manual(values=scen_color_def, labels = function(x) str_wrap(x, width = 24))+
    ggpubr::theme_pubr()
  
  cost_df <- costs %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), scen_list)) %>%
    mutate(name_num = as.numeric(model)) %>%
    filter(model %in% scen_list_i)
  
  cost_plot <- cost_df %>%
    ggplot(aes(x = name_num)) +
    #geom_vline(xintercept = median(unlist(flexible_costs)), linetype = 2) +
    #geom_vline(xintercept = median(unlist(flexible_learning_costs)), linetype = 2) +
    geom_segment(aes(y = lb, yend = ub, x = name_num, xend = name_num, color = model), size = 1) 
  
  if (scen_list[3] %in% scen_list_i) {
    cost_plot <- cost_plot +
      annotate('segment', y = median(unlist(robust_costs)), yend = median(unlist(robust_costs)), x = 2, xend = 6.5, color = 'gray20', linetype = 1) +
      annotate('segment', y = median(unlist(flexible_costs)), yend = median(unlist(flexible_costs)),  x= 3, xend = 6.5, color = 'gray20', linetype = 2) +
      annotate('segment', y = median(unlist(robust_costs)), yend = median(unlist(flexible_costs)), x = 5, xend = 5, arrow = arrow(length = unit(0.2, "cm")))
  }
  
  if (scen_list[4] %in% scen_list_i) {
    cost_plot <- cost_plot +
      annotate('segment', y = median(unlist(flexible_learning_costs)), yend = median(unlist(flexible_learning_costs)), x= 4, xend = 6.5, color = 'gray20', linetype = 2) +
      annotate('segment', y = median(unlist(robust_costs)), yend = median(unlist(flexible_learning_costs)), x = 6, xend = 6, arrow = arrow(length = unit(0.2, "cm")))
  }
  
  cost_plot <- cost_plot +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = median, color = model), size = 3) +
    scale_color_manual(values=scen_color_def, labels = function(x) str_wrap(x, width = 24))+
    #geom_text(aes(label = model, y = ub + 10e6, x = name_num, color = model)) +
    ggpubr::theme_pubr() +
    scale_y_continuous("Cost", labels = scales::unit_format(prefix = "A$", suffix = "M",scale = 1e-6)) +
    coord_cartesian(ylim = c(0, 2.5e8), xlim = c(0,6.5), expand = F) +
    guides(color = 'none')
  cost_plot_vertical <- cost_plot +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  decisions <- list(baseline = cbind(baseline_decisions), 
                    robust = cbind(robust_decisions), 
                    flexible = cbind(flexible_decisions), 
                    flexible_learning = cbind(flexible_learning_decisions)) %>%
    bind_rows(.id = 'model')
  
  cost_plot_horizontal <- cost_plot + 
    coord_flip(ylim = c(0, 2.4e8), xlim = c(6.5,0), expand = F) + 
    scale_y_continuous("Cost", labels = scales::unit_format(prefix = "A$", suffix = "M",scale = 1e-6), 
                       breaks = (0.5:5.5)*0.5e8) +
    scale_x_reverse() +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank())
  
  if (scen_list[3] %in% scen_list_i) {
    cost_plot_vertical <- cost_plot_vertical +
      annotate('text', y = mean(c(median(unlist(flexible_costs)),median(unlist(robust_costs)))), x = 4.9, label = "(i)", hjust = 1)
    
    cost_plot_horizontal <- cost_plot_horizontal +
      annotate('text', y = mean(c(median(unlist(flexible_costs)),median(unlist(robust_costs)))),
               x = 4.9, label = "(i)", vjust = 0)
  }
  
  if (scen_list[4] %in% scen_list_i) {
    cost_plot_vertical <- cost_plot_vertical +
      annotate('text', y = median(unlist(flexible_learning_costs)) - 1e7, x = 5.9, label = "(ii)", hjust = 1)
    cost_plot_horizontal <- cost_plot_horizontal +
      annotate('text', y = mean(c(median(unlist(flexible_learning_costs)),median(unlist(robust_costs)))),
               x = 5.9, label = "(ii)", vjust = 0)
  }
  
  fcn_avg_decisions <- function(decisions, area, summarize = T) {
    id <- decisions$NewPropID
    x_str <- paste0("x", 1:rr)
    y_str <- paste0("sum_y", 1:rr)
    if (summarize) {
      x <- apply(decisions[x_str] * area, 1, mean)
      y <- apply(decisions[y_str] * area / 12, 1, mean)
      data.frame(NewPropID=id,x,y)
    } else {
      x <- decisions[x_str] * area
      y <- decisions[y_str] * area
      cbind(id, x, y)
    }
  }
  
  prop_decisions <- prop_df %>%
    right_join(decisions, by="NewPropID", multiple = 'all') %>%
    pivot_longer(c('x','y'), names_to = 'stage', values_to = 'decision') %>%
    mutate(area = ifelse(stage == 'y' & model == 'flexible_learning', UpperProp * AREA, UpperProp * AREA * decision)) %>%
    mutate(probability = decision) %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), scen_list)) %>%
    filter(decision > 0) %>%
    mutate(stage = factor(stage, c('x','y'), c('Stage 1', 'Stage 2'))) %>%
    filter(model %in% scen_list_i)
  prop_decisions_plot <- prop_decisions %>%
    filter(area > 20) %>%
    ggplot() +
    geom_sf(data = st_transform(nsw_lga_union, st_crs(prop_centroid))) +
    geom_point(aes(x = X, y = Y, color = model, alpha = decision)) +
    scale_color_colorblind7() +
    scale_alpha_continuous()+
    scale_size_continuous("Offer (ha)", range = c(1,4), limits = c(0, NA), breaks = c(1, 1000, 2000, 3000))+
    facet_grid(stage~model, switch = 'y', labeller = label_wrap_gen()) +
    guides(alpha = 'none', fill = 'none', color = 'none') +
    theme_void() +
    theme(strip.text = element_text(margin = margin(3,3,3,3, "pt")))
  
  avg_decisions <- list(baseline = fcn_avg_decisions(baseline_decisions, baseline_area), 
                        robust = fcn_avg_decisions(robust_decisions, baseline_area), 
                        flexible = fcn_avg_decisions(flexible_decisions, baseline_area), 
                        flexible_learning = fcn_avg_decisions(flexible_learning_decisions, baseline_area)) %>%
    bind_rows(.id = 'model')
  
  avg_decisions_full <- list(baseline = fcn_avg_decisions(baseline_decisions, baseline_area, summarize = F), 
                             robust = fcn_avg_decisions(robust_decisions, baseline_area, summarize = F), 
                             flexible = fcn_avg_decisions(flexible_decisions, baseline_area, summarize = F), 
                             flexible_learning = fcn_avg_decisions(flexible_learning_decisions, baseline_area, summarize = F)) %>%
    bind_rows(.id = 'model')
  colnames(avg_decisions_full)[2] <- "NewPropID"
  
  avg_decisions_lga <- avg_decisions %>%
    left_join(prop_lga_lookup, by = 'NewPropID') %>%
    filter(!is.na(NSW_LGA__2)) %>%
    ungroup() %>%
    group_by(model, NSW_LGA__2) %>%
    summarize(stage1 = sum(x), stage2 = sum(y)) %>%
    pivot_longer(c('stage1', 'stage2'), names_to = 'stage', values_to = 'area') %>%
    filter(area > 0) %>%
    group_by(model, stage) %>%
    mutate(proportion = area / sum(area)) %>%
    mutate(stage = factor(stage, c('stage1', 'stage2'), c('Stage 1', 'Stage 2'))) %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), scen_list))
  
  avg_decisions_lga_model <- nsw_lga %>%
    left_join(avg_decisions_lga, by = 'NSW_LGA__2') %>%
    filter(!is.na(model))
  
  prop_cloropleth <- avg_decisions_lga_model %>%
    ggplot() +
    geom_sf(data=nsw_lga, lwd = 0.2, color = 'gray50') +
    geom_sf(aes(fill = proportion), lwd = 0.2, color = 'gray50') +
    facet_grid(stage~model, switch = 'y', labeller = label_wrap_gen()) +
    scale_fill_viridis_c("Proportion") +
    theme_void()
  
  ## Euler sets
  baseline_robust <- fcn_decision_set(baseline_decisions, robust_decisions, baseline_area) %>% colMeans()
  robust_flexible <- fcn_decision_set(robust_decisions, flexible_decisions, baseline_area) %>% colMeans()
  robust_flexible_learning <- fcn_decision_set(robust_decisions, flexible_learning_decisions, baseline_area) %>% colMeans()
  
  fcn_plot_euler <- function(set,labs=c("A","B"),fills=c("gray80", "gray80")) {
    a <- round(as.numeric(set['a']))
    b <- round(as.numeric(set['b']))
    u <- round(as.numeric(set['u']))
    fit <- euler(c("A" = a, "B" = b, "A&B" = u))
    plot(fit, fills = fills, 
         labels = NULL, 
         quantities = list(type='percent', alpha = c(0,0,1)),
         #legend = list(labels = labs, position = 'bottom')
         )    
  }
  
  euler1 <- fcn_plot_euler(baseline_robust, labs=scen_list[c(1,2)], fills = colorpal[c(1,2)])
  euler2 <- fcn_plot_euler(robust_flexible, labs=scen_list[c(2,3)], fills = colorpal[c(2,3)])
  euler3 <- fcn_plot_euler(robust_flexible_learning, labs=c(scen_list[c(2)], "Flexible &\nLearning"), fills = colorpal[c(2,4)])
  
  maps_plot <-  aus_plot + prop_decisions_plot + plot_layout(widths = c(1,2))
  #ggsave("plots/map_plot.png", maps_plot, width = 3000, height = 1300, units = 'px')
  
  layout <- "
  ABB
  CDE
  CFG
  "
  
  plot1 <- wrap_elements(full=aus_plot) + wrap_elements(full=prop_decisions_plot) + cost_plot_vertical + covenanted_area_plot + covenanted_area_end_plot + year_trend_plot + end_range_plot + 
    plot_layout(design = layout, heights = c(2,1,1.5), widths = c(1,1.2,0.4)) & plot_annotation(tag_levels = 'a') 
  
  layout_1a <- "
  ABB
  CCC
  "
  
  eulers <- as.ggplot(euler1) + as.ggplot(euler2) + as.ggplot(euler3) + cowplot::get_legend(end_range_plot) + plot_layout(nrow = 1)
  plot1a <- wrap_elements(full=aus_plot) + wrap_elements(full=prop_decisions_plot) + wrap_elements(eulers) +
    plot_layout(design = layout_1a, heights = c(1,0.5)) & plot_annotation(tag_levels = 'a')
  
  layout_1b <- "
  AB#
  ACD
  "
  
  plot1b <- cost_plot_vertical + covenanted_area_stages_plot + year_trend_plot + end_range_plot + 
    plot_layout(design = layout_1b, heights = c(1,1), widths = c(1,1.2,0.4), guides = 'collect') & plot_annotation(tag_levels = list(c('a','b','c',''))) & theme(legend.position='bottom')
  
  layout_1b2 <- "
  AB
  CD
  "
  
  plot1b2 <- year_split_inset + end_range_plot + cost_plot_horizontal + guide_area()+
    plot_layout(design = layout_1b2, heights = c(1,1), widths = c(1,0.2), guides = 'collect') & plot_annotation(tag_levels = list(c('a','b','c','')))
  
  # Save plots to list
  plot_list[[i]] <- list(
    aus_plot = aus_plot,
    prop_decisions_plot = prop_decisions_plot,
    cost_plot = cost_plot,
    year_trend_plot = year_trend_plot,
    covenanted_area_plot = covenanted_area_plot,
    covenanted_area_end_plot = covenanted_area_end_plot,
    end_range_plot = end_range_plot,
    plot1 = plot1,
    plot1a = plot1a,
    plot1b = plot1b,
    plot1b2 = plot1b2
  )
}

# Save plots ------

ggsave("plots/plot1.png", plot_list[[4]]$plot1, width = 2800, height = 2300, units = 'px')
ggsave("plots/plot1a.png", plot_list[[4]]$plot1a, width = 2800, height = 1800, units = 'px')
ggsave("plots/plot1b.png", plot_list[[4]]$plot1b, width = 2800, height = 1500, units = 'px')
ggsave("plots/plot1b2.png", plot_list[[4]]$plot1b2, width = 3000, height = 1800, units = 'px')

ggsave("plots/plot1.pdf", plot_list[[4]]$plot1, width = 2800, height = 2300, units = 'px')
ggsave("plots/plot1a.pdf", plot_list[[4]]$plot1a, width = 2800, height = 1800, units = 'px')
ggsave("plots/plot1b.pdf", plot_list[[4]]$plot1b, width = 2800, height = 1500, units = 'px')
ggsave("plots/plot1b2.pdf", plot_list[[4]]$plot1b2, width = 3000, height = 1800, units = 'px')
saveRDS(plot_list, file = "plots/plot_list.rds")

## Flexibility types plot --------
recourse_types <- c('(A)', '(E)', '(A/E)')
flexibility_diff_ns <- flexibility_differences(loop_vec = c(1,3,12), dep_var = 5) 
annotation1 <- c("No Learning", 1, -filter(flexibility_diff_ns, recourse == 'ar'& name == 12)$ub - 0.02, "V(F)")
annotation2 <- c("Full Learning", 1, -filter(flexibility_diff_ns, recourse == 'ar'& name == 1)$ub - 0.02, "V(F+L)")
annotations <- as.data.frame(rbind(annotation1, annotation2))
names(annotations) <- c('name', 'model', 'y', 'label')
annotations$model <- as.numeric(annotations$model)
annotations$y <- as.numeric(annotations$y)
annotations$name <- factor(annotations$name, c('No Learning', 'Partial Learning', 'Full Learning'), c('No Learning', 'Partial Learning', 'Full Learning'))

bar_width = 0.2
plot_list_ns <- list()
for (i in 1:length(recourse_types)) {
  ns_plot <- flexibility_diff_ns %>%
    mutate(median = -median, lb = -lb, ub = -ub) %>% # Plot change in positive axis
    mutate(recourse = factor(recourse, c('ar', 'tr', 'fr'), c('(A)', '(E)', '(A/E)'))) %>%
    mutate(name = factor(name, c(12, 3, 1), c('No Learning', 'Partial Learning', 'Full Learning'))) %>%
    mutate(model = as.numeric(recourse)) %>%
    filter(recourse %in% recourse_types[1:i]) %>%
    ggplot(aes( x = model)) +
    geom_hline(yintercept =  0, color = 'gray70') +
    geom_segment(aes(color = recourse, y = median, yend = median, x = model - bar_width, xend = model + bar_width), size = 1) +
    geom_rect(aes(fill = recourse, ymin = lb, ymax = ub, xmin = model - bar_width, xmax = model + bar_width), alpha = 0.5) +
    geom_text(aes(color = recourse, label = recourse, y = lb + 0.04)) +
    geom_text(data = annotations, aes(x = model, y = y, label = label), size = 3) +
    facet_grid(cols = vars(name), labeller = label_wrap_gen()) +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    ggpubr::theme_pubr() +
    guides(color = 'none', fill = 'none') +
    scale_y_continuous("Cost reduction \n(relative to Inflexible)", labels = scales::unit_format(scale = 100,unit = "%")) +
    coord_cartesian(xlim = c(.8,3.2)) +
    theme(axis.title.x = element_blank(),
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(fill = NA, linewidth = 1))
  
  plot_list_ns[[i]] <- ns_plot
}

saveRDS(plot_list_ns, file = "plots/plot_list_ns.rds")
ggsave("plots/plot2.png", plot_list_ns[[3]], width = 2000, height = 800, units = 'px')

# Sensitivity plots ------

plot_shaded_diff <- function(summary_diff_pct) {
  summary_diff_pct %>%
    as.data.frame() %>%
    mutate(name = as.numeric(name)) %>%
    mutate(recourse = factor(recourse, c('ar', 'tr', 'fr'), c('(A)', '(E)', '(A/E)'))) %>%
    mutate(median = -median, lb = -lb, ub = -ub) %>%
    ggplot(aes(x = name, y = median)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = recourse), alpha = 0.25) +
    geom_hline(yintercept = 0) +
    geom_line(aes(color = recourse)) +
    #scale_y_continuous("Change in conservation cost", labels = scales::unit_format(scale = 100,unit = "%")) +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    coord_cartesian(expand = T) +
    ggpubr::theme_pubr() +
    theme(legend.title = element_blank())
}

plot_line_diff <- function(summary_diff_pct, dodge_width = 0.0005, connect = T) {
  p <- summary_diff_pct %>%
    as.data.frame() %>%
    mutate(name = as.numeric(name)) %>%
    mutate(recourse = factor(recourse, c('ar', 'tr', 'fr'), c('(A)', '(E)', '(A/E)'))) %>%
    mutate(median = -median, lb = -lb, ub = -ub) %>%
    mutate(dodge_pos = name+ (as.numeric(recourse)-2)*dodge_width) %>%
    ggplot(aes(x = dodge_pos, y = median)) +
    geom_errorbar(aes(ymin = lb, ymax = ub, color = recourse), width = dodge_width*2) +
    geom_hline(yintercept = 0) 
  
  if (connect) {
    p <- p + geom_line(aes(color = recourse))
  }

  p +
    geom_point(aes(color = recourse)) +
    scale_y_continuous("Cost reduction", labels = scales::unit_format(scale = 100,unit = "%")) +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    ggpubr::theme_pubr() +
    theme(legend.title = element_blank())
}

plot_line_diff_learning <- function(summary_diff_pct, dodge_width = 0.0005, label = "", connect = T) {
  p <- summary_diff_pct %>%
    as.data.frame() %>%
    mutate(name = as.numeric(name)) %>%
    mutate(learning = factor(learning, c('no', 'full'), scen_list[3:4])) %>%
    mutate(median = -median, lb = -lb, ub = -ub) %>%
    ggplot(aes(x = name, y = median, color = learning)) +
    
    geom_hline(yintercept = 0) 
  
  if (connect) {
    p <- p + geom_line(position=position_dodge(width=dodge_width))
  }
  
  p +
    geom_point(position=position_dodge(width=dodge_width)) +
    geom_errorbar(aes(ymin = lb, ymax = ub), width = dodge_width, position=position_dodge(width=dodge_width)) + 
    scale_y_continuous("Cost reduction", labels = scales::unit_format(scale = 100,unit = "%")) +
    scale_x_continuous(label) +
    scale_color_manual(values = scen_color_def[3:4]) +
    scale_fill_manual(values = scen_color_def[3:4]) +
    ggpubr::theme_pubr() +
    theme(legend.title = element_blank())
}

# Extract matrices of the differences between robust and flexible solutions -----
full_learning <- c(7000,1,6,0.25,1,30,0.02,0.1)
no_learning <- c(7000,1,6,0.25,12,30,0.02,0.1)

cost_first_stage <- read_csv(paste0(results_dir, 'cost1_', get_run_string(no_learning), '.csv'))

dr_vec <- seq(0.0,1,0.05) %>%
  sapply(function(x) format(round(x, 2), nsmall = 1))
dr_flex_diff <- list(
  no = flexibility_differences(params = no_learning, loop_vec = dr_vec, dep_var = 8),
  full = flexibility_differences(params = full_learning, loop_vec = dr_vec, dep_var = 8)
  ) %>%
  bind_rows(.id = "learning")

dr_flex_plot <- dr_flex_diff %>%
  mutate(learning = factor(learning, c('no', 'full'), scen_list[3:4])) %>%
  ggplot(aes(x = as.numeric(name))) +
  geom_vline(xintercept = 7.2/101, linetype = 'longdash', color = 'gray50') +
  annotate("text", x = .02+7.2/101, y = .95, label = "", hjust = 0) +
  geom_ribbon(aes(fill = learning, ymin = -lb, ymax = -ub), alpha = .4) +
  geom_point(aes(color = learning, y = -median)) +
  geom_line(aes(color = learning, y = -median)) +
  scale_y_continuous("Cost reduction", labels = scales::unit_format(scale = 100,unit = "%"), limits = c(0,1)) +
  scale_color_manual("", values = scen_color_def, labels = function(x) str_wrap(x, width = 24)) +
  scale_fill_manual("", values = scen_color_def, labels = function(x) str_wrap(x, width = 24)) +
  scale_x_continuous("Probability of land clearing") +
  ggpubr::theme_pubr() +
  labs(caption = "Dashed line show 1972-2014 total deforestation in Australia (Evans, 2016)") +
  theme(legend.position = "bottom")

## Plot change in protected area size in stage 1 relative to stage 2
dr_share_full_learning <- lapply(dr_vec, function(i) {
  l <- full_learning
  l[8] <- i
  decision <- read_csv(paste0(results_dir, 'decision_ar_', get_run_string(l), '.csv'), col_types = cols())
  baseline_area <- read_csv(paste0(results_dir, "area_", get_run_string(l), ".csv"), col_types = cols()) 
  x <- colSums(decision[,paste0('x', 1:rr)] * baseline_area)
  y <- colSums(decision[,paste0('sum_y', 1:rr)]/12 * baseline_area)
  x / (x+y)
})
dr_share_no_learning <- lapply(dr_vec, function(i) {
  l <- no_learning
  l[8] <- i
  decision <- read_csv(paste0(results_dir, 'decision_ar_', get_run_string(l), '.csv'), col_types = cols())
  baseline_area <- read_csv(paste0(results_dir, "area_", get_run_string(l), ".csv"), col_types = cols()) 
  x <- colSums(decision[,paste0('x', 1:rr)] * baseline_area)
  y <- colSums(decision[,paste0('sum_y', 1:rr)]/12 * baseline_area)
  x / (x+y)
})

dr_share_func <- function(df_list) {
  a <- sapply(df_list, function(df) {
    v <- unlist(df)
    data.frame(median = median(v), lb = min(v), ub = max(v))
  }) %>%
    t() %>%
    as.data.frame()
  a$dr <- as.numeric(dr_vec)
  a
}

dr_share_diff <- list(
  full = dr_share_func(dr_share_full_learning),
  no = dr_share_func(dr_share_no_learning)
) %>%
  bind_rows(.id = 'learning')

dr_share_plot <- dr_share_diff %>%
  as.data.frame() %>%
  mutate(learning = factor(learning, c('no', 'full'), scen_list[3:4])) %>%
  mutate(median = as.numeric(median), lb = as.numeric(lb), ub = as.numeric(ub)) %>%
  ggplot(aes(color = learning, fill = learning, x = dr, y = median)) +
  #geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.5) +
  geom_point() +
  geom_line() +
  scale_color_manual("", values = scen_color_def) +
  scale_fill_manual("", values = scen_color_def) +
  geom_vline(xintercept = 7.2/101, linetype = 'longdash', color = 'gray50') +
  scale_y_continuous("Median proportion of land \nprotected in Stage 1", labels = scales::unit_format(scale = 100,unit = "%"), limits = c(0,1)) +
  ggpubr::theme_pubr() +
  guides(color = 'none', fill = 'none')+
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

dr_plots <- dr_share_plot / dr_flex_plot + theme(legend.position = 'bottom') & plot_annotation(tag_levels = 'a')

ggsave("plots/plot3.png", dr_plots, width = 1500, height = 2200, units = 'px')

dr_flex_diff %>%
  ggplot(aes(color = learning, fill = learning, x = as.numeric(name))) +
  geom_point(aes(y = -median), position=position_dodge(width=0.02)) +
  geom_errorbar(aes(x=as.numeric(name), ymin = -lb, ymax = -ub), width = 0.02, position=position_dodge(width=0.02)) + 
  scale_y_continuous("Cost reduction", labels = scales::unit_format(scale = 100,unit = "%")) +
  ggsci::scale_color_d3() +
  ggpubr::theme_pubr()

sp_flex_diff <- list(
  full = flexibility_differences(params = full_learning, loop_vec = 1:10, dep_var = 2),
  no = flexibility_differences(params = no_learning, loop_vec = 1:10, dep_var = 2)
) %>%
  bind_rows(.id = 'learning')

tt_flex_diff <- list(
  full = flexibility_differences(params = full_learning, loop_vec = c(3,4,5,6), dep_var = 3, realisation = 1:10),
  no = flexibility_differences(params = no_learning, loop_vec = c(3,4,5,6), dep_var = 3, realisation = 1:10)
) %>%
  bind_rows(.id = 'learning')

kt_flex_diff <- list(
  full = flexibility_differences(params = full_learning, loop_vec = c(0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275), dep_var = 4),
  no = flexibility_differences(params = no_learning, loop_vec = c(0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275), dep_var = 4)
) %>%
  bind_rows(.id = 'learning')

# Social discount rates
sdr_flex_diff <- list(
  no = flexibility_differences(params = no_learning, loop_vec = c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05), dep_var = 7),
  full = flexibility_differences(params = full_learning, loop_vec = c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05), dep_var = 7)
) %>%
  bind_rows(.id = 'learning')

# Policy targets
k_flex_diff <- list(
  full = flexibility_differences(params = full_learning, loop_vec = 1000+1500*(0:7), dep_var = 1),
  no = flexibility_differences(params = no_learning, loop_vec = 1000+1500*(0:7), dep_var = 1)
) %>%
  bind_rows(.id = 'learning')

tt_flex_plot <- plot_line_diff_learning(mutate(tt_flex_diff, name = 2000+as.numeric(name)*10), 0.2, "Year for new investment (t')")
kt_flex_plot <- plot_line_diff_learning(kt_flex_diff, 0.01, "Ecological indicator cut-off")
sdr_flex_plot <- plot_line_diff_learning(sdr_flex_diff, 0.001, "Discount rate (1-ρ)")
k_flex_plot <- plot_line_diff_learning(k_flex_diff, 500, "Policy target (hectares of koala habitat)")

sensitivity_plots <- (k_flex_plot + sdr_flex_plot + tt_flex_plot + kt_flex_plot ) + plot_layout(guides='collect') & theme(legend.position = 'bottom') & plot_annotation(tag_levels = 'a')
ggsave("plots/sensitivity_plots.png", sensitivity_plots, units = 'px', width = 3000, height = 2000)



