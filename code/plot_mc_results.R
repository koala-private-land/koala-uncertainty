## Plots the results from sensitivity analyses

library(tidyverse)
library(ggplot2)
library(purrr)
library(patchwork)
library(geodata)
library(sf)
library(nngeo)

results_dir <- "./results/mc_sim/"
#results_dir <- "~/Documents/OneDrive - The University of Queensland/Documents/GitHub/koala-uncertainty/results/mc_sim/"

flexibility_differences <- function(recourse = NULL, loop_vec = 1:10, params = c(sp, tt, kt, ns, rr), dep_var = 1, interval = c(0.05, 0.95), realisation = NULL, base_dir = results_dir) {
  base_string <- "run_sp-%s_tt-%s_kt-%s_ns-%s_r-%s"
  cost_diff_pct <- lapply(loop_vec, function(x) {
    vec_var <- params
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
  
  return(summary_diff_pct)
}


year_vec <- seq(2000, 2070, by=10)

# Baseline vs robust ------

# Load property data
properties <- st_read("data/spatial_predictions_v1_1.gdb", layer = "spatial_pred_inperp_v1_1")
prop_centroid <- properties %>%
  st_centroid() %>%
  select(NewPropID, MeanProp, UpperProp, Shape_Area) %>%
  mutate(AREA = Shape_Area)
prop_df <- st_drop_geometry(prop_centroid) %>%
  cbind(st_coordinates(prop_centroid))
## Spatial plots
aus_border <- gadm(country = "AUS", level = 1, path = "data/", resolution = 2) %>%
  sf::st_as_sf()
nsw_lga <- st_read( "data/planning_units.gdb", layer = 'nsw_lga_pu_all') %>%
  st_transform(4326) %>%
  st_simplify(preserveTopology = T, dTolerance = 2000)
nsw_lga_union <- nsw_lga%>%
  st_union() %>%
  nngeo::st_remove_holes()
nsw_bbox <- st_bbox(nsw_lga_union)
bbox_buffer <- 1
aus_xlim <- c(113.338953078, 153.569469029)
aus_ylim <- c(-43.6345972634,-10.6681857235)
aus_plot <- ggplot() +
  geom_sf(data = aus_border, fill = 'gray80', color = 'white') +
  geom_rect(aes(xmin = nsw_bbox$xmin - bbox_buffer, xmax = nsw_bbox$xmax + bbox_buffer, ymin = nsw_bbox$ymin - bbox_buffer, ymax = nsw_bbox$ymax + bbox_buffer), fill = NA, color = 'gray50') +
  geom_sf(data = nsw_lga_union, linewidth = 0.5, color = NA, fill = '#FFD580') +
  coord_sf(xlim = aus_xlim, ylim = aus_ylim) +
  theme_void()
#theme_bw() +
#theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
#      panel.background = element_rect(fill = 'aliceblue'))

# Default parameters
sp = 1 #1:10
tt = 6 # [0,1,2,3,4,5,6,7]
kt = 0.25 # [0.1, 0.15, 0.2, 0.25, 0.3]
ns = 12 # 1:12
rr = 10

scen_list <- c('CC', 'RI', 'F', 'F+L')
plot_list <- list()

for (i in 1:length(scen_list)) {
  scen_list_i <- scen_list[1:i]

  get_run_string <- function(param=default) do.call(sprintf, c(fmt = "run_sp-%s_tt-%s_kt-%s_ns-%s_r-%s", as.list(param)))
  baseline_string <- get_run_string(c(sp, tt, kt, ns, rr))
  robust_string   <- get_run_string(c(sp, tt, kt, ns, rr))
  flexible_string <- get_run_string(c(sp, tt, kt, ns, rr))
  flexible_learning_string <- get_run_string(c(sp, tt, kt, 1, rr))
  
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
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), c('CC', 'RI', 'F', 'F+L'))) 
  
  year_trend_plot <- baseline_robust_df %>%
    filter(year >= 2020) %>%
    filter(model %in% scen_list_i) %>%
    ggplot(aes(x=year)) +
    geom_hline(yintercept = 7000) +
    geom_vline(xintercept = 2020, linetype = 2, color = 'gray50') +
    geom_vline(xintercept = year_vec[tt], linetype = 2, color = 'gray50') +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = model), alpha = 0.25) +
    geom_line(aes(y = median, color = model), linewidth = 1) +
    geom_point(data = filter(baseline_robust_df, t %in% c(3)), aes(y = median, color = model), size = 2)+
    geom_point(data = filter(baseline_robust_df, t %in% tt & model %in% c('F','F+L')), aes(y = median, color = model), size = 2)+
    scale_y_continuous("Koala habitat (ha)") +
    scale_x_continuous("Year") +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    coord_cartesian(xlim = c(2020,2070), ylim = c(0, 15000)) +
    guides(color = 'none', fill = 'none')+
    ggpubr::theme_pubr()
  
  end_range_plot <- baseline_robust_df %>%
    filter(year == 2070) %>%
    mutate(name_num = as.numeric(as.factor(model))) %>%
    filter(model %in% scen_list_i) %>%
    ggplot(aes(x = name_num, y = median, fill = model)) +
    geom_hline(yintercept = 7000) +
    geom_rect(aes(ymin = lb, ymax = ub, xmin = name_num - 0.4, xmax = name_num + 0.4))+
    geom_segment(aes(y = median, yend = median, x = name_num-0.4, xend = name_num+0.4), color = 'white', linewidth = 1) +
    geom_text(aes(y = 0, x = name_num, label = model, color = model), size = 4, vjust = 0) +
    coord_cartesian(ylim = c(0, 15000)) +
    guides(color = 'none', fill = 'none') +
    ggsci::scale_fill_nejm() +
    ggsci::scale_color_nejm() +
    theme_void() +
    theme(axis.line.x = element_blank())
  
  ## Import decision vectors
  baseline_decisions <- read_csv(paste0(results_dir, "decision_baseline_", baseline_string, ".csv"), col_types = cols())  %>% mutate(y = sum_y / 12, w = sum_w / 12)
  robust_decisions <- read_csv(paste0(results_dir, "decision_nr_", robust_string, ".csv"), col_types = cols()) %>% mutate(y = sum_y / 12, w = sum_w / 12)
  flexible_decisions <- read_csv(paste0(results_dir, 'decision_ar_', flexible_string, ".csv"), col_types = cols()) %>% mutate(y = sum_y / 12, w = sum_w / 12)
  flexible_learning_decisions <- read_csv(paste0(results_dir, 'decision_ar_', flexible_learning_string, ".csv"), col_types = cols()) %>% mutate(y = sum_y / 12, w = sum_w / 12)
  
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
  
  ## Sum of area
  spatial_pred <- read_csv('data/spatial_predictions_10yr.csv')
  calculate_offered_area <- function(spatial_pred, decisions) {
    spatial_pred_subset <- right_join(spatial_pred, decisions, by = 'NewPropID')
    expected_area <- spatial_pred_subset %>%
      mutate(area = UpperProp * AREA) %>%
      select(area)
    return(as.data.frame(t(expected_area$area) %*% as.matrix(decisions[,c('x','y','w')])))
  }
  calculate_area <- function(area, decisions) {
    t(as.matrix(area)) %*% as.matrix(decisions)
  }
  
  area_range <- function(area) {
    area <- as.data.frame(area)
    median = apply(area, 2, quantile, probs = .5)
    lb = area[which.min(area$x + area$y - area$w),]
    ub = area[which.max(area$x + area$y - area$w),]
    df <- rbind(median = median, lb = lb, ub = ub) %>% as.data.frame()
    df$stat <- as.list(row.names(df))
    return(df)
  }
  
  baseline_sum_area <- calculate_area(baseline_area, baseline_decisions[,c('x','y','w')]) %>% area_range()
  robust_sum_area <- calculate_area(robust_area, robust_decisions[,c('x','y','w')]) %>% area_range()
  flexible_sum_area <- calculate_area(flexible_area, flexible_decisions[,c('x','y','w')]) %>% area_range()
  flexible_learning_sum_area <- calculate_area(flexible_learning_area, flexible_learning_decisions[,c('x','y','w')]) %>% area_range()
  
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
  
  baseline_offered_area <- calculate_offered_area(spatial_pred, baseline_decisions) %>% area_to_ts()
  robust_offered_area <- calculate_offered_area(spatial_pred, robust_decisions)%>% area_to_ts()
  flexible_offered_area <- calculate_offered_area(spatial_pred, flexible_decisions)%>% area_to_ts()
  flexible_learning_offered_area <- calculate_offered_area(spatial_pred, flexible_learning_decisions)%>% area_to_ts()
  
  offered_area <- list(baseline = baseline_offered_area, robust = robust_offered_area, 
                       flexible = flexible_offered_area, flexible_learning = flexible_learning_offered_area) %>%
    bind_rows(.id = 'model') %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), c('CC', 'RI', 'F', 'F+L'))) %>%
    filter(year >= 2020)
  
  covenanted_area_plot <- list(baseline = baseline_area_ts, robust = robust_area_ts, flexible = flexible_area_ts, flexible_learning = flexible_learning_area_ts) %>%
    bind_rows(.id = 'model') %>%
    pivot_wider(names_from = 'stat', values_from = 'area') %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), c('CC', 'RI', 'F', 'F+L'))) %>%
    filter(year >= 2020) %>%
    filter(model %in% scen_list_i) %>%
    ggplot(aes(x = year)) +
    geom_vline(xintercept = 2020, linetype = 2, color = 'gray50') +
    geom_vline(xintercept = year_vec[tt], linetype = 2, color = 'gray50') +
    geom_ribbon(aes(fill = model, ymin = lb, ymax = ub), alpha = 0.25) +
    geom_line(data = offered_area, aes(color = model, y = area), linewidth = 1) +
    #geom_line(aes(color = model, y = median)) +
    
    scale_y_continuous("Covenant offers/ \nyield (ha)") +
    scale_x_continuous("Year") +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    coord_cartesian(xlim = c(2020,2070), ylim = c(4000, 18000)) +
    guides(color = 'none', fill = 'none')+
    ggpubr::theme_pubr() +
    annotate('text', x = mean(c(year_vec[3], year_vec[tt])), y = 18000, color = 'black', label = 'Stage 1') +
    annotate('text', x = mean(c(year_vec[tt], year_vec[length(year_vec)])) , y = 18000, color = 'black', label = 'Stage 2') +
    theme(axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  offered_area_end <- offered_area %>%
    mutate(name_num = as.numeric(model)) %>%
    filter(year == 2070)
  
  covenanted_area_end_plot <- list(baseline = baseline_area_ts, robust = robust_area_ts, flexible = flexible_area_ts, flexible_learning = flexible_learning_area_ts) %>%
    bind_rows(.id = 'model') %>%
    pivot_wider(names_from = 'stat', values_from = 'area') %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), c('CC', 'RI', 'F', 'F+L'))) %>%
    filter(year == 2070) %>%
    filter(model %in% scen_list_i) %>%
    mutate(name_num = as.numeric(as.factor(model))) %>%
    ggplot(aes(x = name_num, y = median, fill = model)) +
    #geom_hline(yintercept = 7000) +
    geom_rect(aes(ymin = lb, ymax = ub, xmin = name_num - 0.4, xmax = name_num + 0.4))+
    geom_point(data = offered_area_end, aes(y = area, x = name_num, color = model), size = 2, shape = 23) +
    geom_text(aes(y = 5000, x = name_num, label = model, color = model), size = 4, vjust = 0) +
    coord_cartesian(ylim = c(4000, 18000)) +
    guides(color = 'none', fill = 'none') +
    ggsci::scale_fill_nejm() +
    ggsci::scale_color_nejm() +
    theme_void() +
    theme(axis.line.x = element_blank())
  
  ## Cost plot ----------
  costs <- list(baseline = baseline_costs, robust = robust_costs, flexible = flexible_costs, flexible_learning = flexible_learning_costs)%>%
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
  
  cost_df <- costs %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), c('CC', 'RI', 'F', 'F+L'))) %>%
    mutate(name_num = as.numeric(model))
  
  cost_plot <- cost_df %>%
    ggplot(aes(x = name_num)) +
    #geom_vline(xintercept = median(unlist(flexible_costs)), linetype = 2) +
    #geom_vline(xintercept = median(unlist(flexible_learning_costs)), linetype = 2) +
    geom_segment(aes(y = lb, yend = ub, x = name_num, xend = name_num, color = model), linewidth = 1) +
    annotate('segment', y = median(unlist(robust_costs)), yend = median(unlist(robust_costs)), x = 2, xend = 6.5, color = 'gray20', linetype = 1) +
    annotate('segment', y = median(unlist(flexible_costs)), yend = median(unlist(flexible_costs)),  x= 3, xend = 6.5, color = 'gray20', linetype = 2) +
    annotate('segment', y = median(unlist(flexible_learning_costs)), yend = median(unlist(flexible_learning_costs)), x= 4, xend = 6.5, color = 'gray20', linetype = 2) +
    annotate('text', y = mean(c(median(unlist(flexible_costs)),median(unlist(robust_costs)))), x = 5.5, label = "V(F)", hjust = 1) +
    annotate('text', y = mean(c(median(unlist(flexible_costs)),median(unlist(flexible_learning_costs)))), x = 5.9, label = "V(F+L)", hjust = 1) +
    annotate('segment', y = median(unlist(robust_costs)), yend = median(unlist(flexible_costs)), x = 5.6, xend = 5.6, arrow = arrow(length = unit(0.2, "cm"))) +
    annotate('segment', y = median(unlist(robust_costs)), yend = median(unlist(flexible_learning_costs)), x = 6, xend = 6, arrow = arrow(length = unit(0.2, "cm"))) +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = median, color = model), size = 3) +
    ggsci::scale_color_nejm()+
    geom_text(aes(label = model, y = lb - 5e6, x = name_num, color = model), size = 6) +
    ggpubr::theme_pubr() +
    scale_y_continuous("Cost", labels = scales::unit_format(prefix = "A$", suffix = "M",scale = 1e-6)) +
    coord_cartesian(ylim = c(0, 1.2e8), xlim = c(0,6.5), expand = F) +
    #scale_y_reverse() +
    guides(color = 'none') +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  decisions <- list(baseline = baseline_decisions, robust = robust_decisions, flexible = flexible_decisions, flexible_learning = flexible_learning_decisions) %>%
    bind_rows(.id = 'model')
  prop_decisions_plot <- prop_df %>%
    right_join(decisions) %>%
    pivot_longer(c('x','y'), names_to = 'stage', values_to = 'decision') %>%
    mutate(area = ifelse(stage == 'y' & model == 'flexible_learning', UpperProp * AREA,UpperProp * AREA * decision)) %>%
    mutate(probability = ifelse(stage == 'y' & model == 'flexible_learning', decision, decision > 0)) %>%
    mutate(model = factor(model, c('baseline', 'robust', 'flexible', 'flexible_learning'), c('CC', 'RI', 'F', 'F+L'))) %>%
    filter(decision > 0) %>%
    mutate(stage = factor(stage, c('x','y'), c('Stage 1', 'Stage 2'))) %>%
    filter(model %in% scen_list_i) %>%
    ggplot() +
    geom_sf(data = st_transform(nsw_lga_union, st_crs(prop_centroid))) +
    geom_point(aes(x = X, y = Y, size = area, fill = model, alpha = probability), color = 'white', pch = 21) +
    ggsci::scale_fill_nejm() +
    facet_grid(stage~model, switch = 'y') +
    guides(size = 'none', alpha = 'none', fill = 'none') +
    theme_void()
  prop_decisions_plot
  
  maps_plot <-  aus_plot + prop_decisions_plot + plot_layout(widths = c(1,2))
  #ggsave("plots/map_plot.png", maps_plot, width = 3000, height = 1300, units = 'px')
  
  layout <- "
  ABB
  CDF
  CEG
  "
  
  plot1 <- wrap_elements(full=aus_plot) + wrap_elements(full=prop_decisions_plot) + cost_plot + covenanted_area_plot + year_trend_plot + covenanted_area_end_plot+ end_range_plot + 
    plot_layout(design = layout, heights = c(2,1,1.5), widths = c(1,1.2,0.4)) & plot_annotation(tag_levels = 'a') 

  ## Flexibility types plot --------
  
  flexibility_diff_ns <- flexibility_differences(loop_vec = c(1,3,12), dep_var = 4) 
  annotation1 <- c("No Learning", 1, -filter(flexibility_diff_ns, recourse == 'ar'& name == 12)$ub - 0.02, "V(F)")
  annotation2 <- c("Full Learning", 1, -filter(flexibility_diff_ns, recourse == 'ar'& name == 1)$ub - 0.02, "V(F+L)")
  annotations <- as.data.frame(rbind(annotation1, annotation2))
  names(annotations) <- c('name', 'model', 'y', 'label')
  annotations$model <- as.numeric(annotations$model)
  annotations$y <- as.numeric(annotations$y)
  annotations$name <- factor(annotations$name, c('No Learning', 'Partial Learning', 'Full Learning'), c('No Learning', 'Partial Learning', 'Full Learning'))
  
  bar_width = 0.2
  ns_plot <- flexibility_diff_ns %>%
    mutate(median = -median, lb = -lb, ub = -ub) %>% # Plot change in positive axis
    mutate(recourse = factor(recourse, c('ar', 'tr', 'fr'), c('(A)', '(E)', '(A/E)'))) %>%
    mutate(name = factor(name, c(12, 3, 1), c('No Learning', 'Partial Learning', 'Full Learning'))) %>%
    mutate(model = as.numeric(recourse)) %>%
    ggplot(aes( x = model)) +
    geom_hline(yintercept =  0, color = 'gray70') +
    geom_segment(aes(color = recourse, y = median, yend = median, x = model - bar_width, xend = model + bar_width), linewidth = 1) +
    geom_rect(aes(fill = recourse, ymin = lb, ymax = ub, xmin = model - bar_width, xmax = model + bar_width), alpha = 0.5) +
    geom_text(aes(color = recourse, label = recourse, y = lb + 0.04)) +
    geom_text(data = annotations, aes(x = model, y = y, label = label), size = 3) +
    facet_grid(cols = vars(name)) +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    ggpubr::theme_pubr() +
    guides(color = 'none', fill = 'none') +
    scale_y_continuous("Cost reduction \n(relative to Robust)", labels = scales::unit_format(scale = 100,unit = "%")) +
    theme(axis.title.x = element_blank(),
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(fill = NA, linewidth = 1))
  
  # Save plots to list
  plot_list[[i]] <- list(
    aus_plot = aus_plot,
    prop_decisions_plot = prop_decisions_plot,
    cost_plot = cost_plot,
    year_trend_plot = year_trend_plot,
    covenanted_area_plot = covenanted_area_plot,
    covenanted_area_end_plot = covenanted_area_end_plot,
    end_range_plot = end_range_plot, 
    ns_plot = ns_plot,
    plot1 = plot1
  )
}

# Save plots ------

ggsave("plots/plot1.png", plot_list[[4]]$plot1, width = 2800, height = 2300, units = 'px')
ggsave("plots/plot2.png", plot_list[[4]]$ns_plot, width = 2000, height = 800, units = 'px')

library(gganimate)
year_trend_anim <- year_trend_plot +
  transition_states(as.numeric(model))+
  shadow_mark(alpha = alpha/2)   
anim_save(year_trend_anim)

# Sensitivity plots ------

plot_shaded_diff <- function(summary_diff_pct) {
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

default <- c(1,5,0.25,12,10)
flexibility_differences(params = default, loop_vec = 1:10, dep_var = 1) %>% plot_shaded_diff()
flexibility_differences(params = default, loop_vec = c(2,3,4,5, 6,7), dep_var = 2, realisation = 1:10)%>% plot_shaded_diff()
flexibility_differences(params = default, loop_vec = c(0.1, 0.15, 0.2, 0.25), dep_var = 3, realisation = 1:10)%>% plot_shaded_diff()
flexibility_differences(params = default, loop_vec = 1:12, dep_var = 4, realisation = 1:10, interval = c(0.05, 0.95)) %>% plot_shaded_diff() +
  scale_x_reverse("Learning (number of climate scenarios)") +
  scale_y_reverse("Cost reduction", labels = scales::unit_format(scale = 100,unit = "%"), limits = c(-1, 0.1))

