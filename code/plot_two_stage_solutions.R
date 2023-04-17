library(sf)
library(tidyverse)
library(ggplot2)

properties_subset <- st_read("data/spatial_predictions_v1_1.gdb/", layer="spatial_pred_10yr_v1_1")
properties_centroid <- st_centroid(properties_subset)

nsw_lga <- st_read( "data/planning_units.gdb", layer = 'nsw_lga_pu_all') %>%
  st_transform(4326) %>%
  st_simplify(preserveTopology = T, dTolerance = 2000)

solution_list <- list()
solution_list$no_recourse    <- read_csv('data/solution_no_recourse_kitl.csv')
solution_list$add_only       <- read_csv('data/solution_add_only_kitl.csv')
solution_list$terminate_only <- read_csv('data/solution_terminate_only_kitl.csv')
solution_list$full_recourse <- read_csv('data/solution_full_recourse_kitl.csv')
factor_recourse_var <- function(r) factor(r, c("no_recourse", "add_only", "terminate_only", "full_recourse"),
       c("No recourse", "Add only", "Terminate only", "Full recourse"))
solution <- solution_list %>% 
  bind_rows(.id = "recourse") %>%
  mutate(recourse = factor_recourse_var(recourse))

properties_solution <- properties_centroid %>%
  right_join(solution)

obj_values <- read_csv("data/obj_values.csv")

obj_values_plot <- obj_values %>%
  mutate(Solution = factor(Solution, c("no_recourse", "add_only", "terminate_only", "full_recourse"),
                           c("No recourse", "Add only", "Terminate only", "Full recourse"))) %>%
  ggplot(aes(y = `Objective Value`, x = fct_rev(Solution))) +
  geom_bar(stat="identity", width = 0.5) +
  geom_text(aes(label = paste(round(`Objective Value`/1e6,1),"M")), nudge_y = -3e6, color = 'white') +
  scale_x_discrete('') +
  scale_y_continuous("Total conservation costs (AUD)", labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  ggpubr::theme_pubclean()
obj_values_plot

first_stage <- ggplot(properties_solution) +
  geom_sf(data = nsw_lga, fill = 'white') +
  geom_sf(data = properties_solution %>% filter(x == 0), color = 'gray60', size = 0.1) +
  geom_sf(data = properties_solution %>% filter(x == 1), color = '#20854EFF') +
  facet_wrap(~recourse, nrow = 1) +
  theme_void()
first_stage

properties_solution_long <- properties_solution %>%
  pivot_longer(cols = c("sum_y", "sum_w"), names_to = "variable", values_to = "decision") %>%
  mutate(variable = factor(variable, c("sum_y", "sum_w"), c("Add", "Terminate")))

second_stage <- ggplot(properties_solution_long) +
  geom_sf(data = nsw_lga, fill = 'white') +
  geom_sf(data = properties_solution_long %>% filter(decision < 1), color = 'gray60', size = 0.1) +
  geom_sf(data = properties_solution_long %>% filter(decision >= 1), aes(color = decision)) +
  scale_color_viridis_b("", direction = -1) +
  facet_grid(rows = vars(variable), cols = vars(recourse)) +
  theme_void()
second_stage

ggsave("plots/first_stage_solution.png", first_stage)
ggsave("plots/second_stage_solution.png", second_stage)

obj_values <- list()
obj_values$no_recourse <- read_csv("data/cost_no_recourse.csv")
obj_values$add_only <- read_csv("data/cost_add_recourse.csv")
obj_values$terminate_only <- read_csv("data/cost_terminate_recourse.csv")
obj_values$full_recourse <- read_csv("data/cost_full_recourse.csv")
obj_values_df <- obj_values %>%
  bind_rows() %>%
  t()
colnames(obj_values_df) <-c('no_recourse', 'add_recourse', 'terminate_recourse', 'full_recourse')
obj_values_plot <- obj_values_df %>%
  as.data.frame() %>%
  pivot_longer(colnames(obj_values_df)) %>%
  mutate(name = factor(name, colnames(obj_values_df), c('None', 'Add', 'Terminate', 'Full'))) %>%
  ggplot(aes(x = value, y = name)) +
  geom_vline(xintercept = obj_values_df[1,1], color = 'gray50') +
  geom_boxplot() +
  scale_x_continuous("Conservation cost", labels = scales::unit_format(prefix = "A$", suffix = "M",scale = 1e-6)) +
  scale_y_discrete('') +
  ggpubr::theme_pubr()
ggsave("plots/obj_values.png", obj_values_plot)

pct_saving_plot <- obj_values_df %>%
  as.data.frame() %>%
  mutate(pct_change = (add_recourse - no_recourse)/no_recourse,
         climate_model = rownames(obj_values_df)) %>%
  mutate(climate_model = factor(climate_model, labels = rownames(obj_values_df), levels = rownames(obj_values_df)[order(pct_change)])) %>%
  ggplot(aes(x = pct_change, y = climate_model, fill = pct_change)) +
  geom_bar(stat = 'identity') +
  scale_y_discrete() +
  scale_x_continuous("Conservation costs", labels = scales::percent, limits = c(-.5,0.1)) +
  #scale_x_continuous(labels = scales::unit_format(suffix = "M",scale = 1e-6)) +
  geom_vline(xintercept = 0) +
  guides(fill = 'none')+
  theme_minimal()
ggsave("plots/pct_savings.png", pct_saving_plot)

metric <- list()
metric$no_recourse <- read_csv("data/habitat_size_no_recourse.csv")
metric$add_only <- read_csv("data/habitat_size_add_recourse.csv")

metric_df <- metric %>%
  lapply(function(x) {
    list(mean = rowMeans(x), lb = apply(x, 1, min), ub = apply(x, 1, max), t = 2000+10*(0:7))
  }) %>%
  bind_rows(.id = 'model')

habitat_size_plot <- metric_df %>%
  mutate(model = factor(model, levels = c('no_recourse','add_only'), labels = c('No recourse', 'Add recourse')))%>%
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = 'gray70') +
  geom_line(aes(y = mean)) +
  geom_hline(yintercept = 7000) +
  scale_y_continuous("Protected koala habitat (ha)", limits = c(6000, 12000)) +
  scale_x_continuous("Year") +
  facet_wrap(~model) +
  ggpubr::theme_pubr()
ggsave("plots/habitat_size_plot.png", habitat_size_plot)

saa_metric <- list()
saa_metric$nr_nsaa <- read_csv("data/nr_metric.csv") %>%
  mutate(t = 2000+10*(0:7))
saa_metric$add_nsaa <- read_csv("data/add_metric.csv")%>%
  mutate(t = 2000+10*(0:7))
saa_metric$nr_saa <- read_csv("data/saa_nr_metric.csv")%>%
  mutate(t = 2000+10*(0:7))
saa_metric$add_saa <- read_csv("data/saa_add_metric.csv")%>%
  mutate(t = 2000+10*(0:7))
saa_p1 <- saa_metric %>%
  bind_rows(.id = 'name') %>%
  separate(name, c('model', 'saa')) %>%
  mutate(model = factor(model, levels = c('nr','add'), labels = c('No recourse', 'Add recourse')))%>%
  filter(saa == 'nsaa') %>%
  mutate(saa = factor(saa, levels = c('saa','nsaa'), labels = c('Adoption uncertainty', 'No adoption uncertainty')))%>%
  ggplot(aes( x = t, fill = saa)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.4) +
  geom_hline(yintercept = 7000) +
  geom_line(aes(y = mean, color = saa)) +
  scale_y_continuous("Protected koala habitat (ha)", limits = c(5000, 12000)) +
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  facet_wrap(~model) +
  ggpubr::theme_pubr()
saa_p2 <- saa_metric %>%
  bind_rows(.id = 'name') %>%
  separate(name, c('model', 'saa')) %>%
  mutate(model = factor(model, levels = c('nr','add'), labels = c('No recourse', 'Add recourse')))%>%
  #filter(saa == 'nsaa') %>%
  mutate(saa = factor(saa, levels = c('saa','nsaa'), labels = c('Adoption uncertainty', 'No adoption uncertainty')))%>%
  ggplot(aes( x = t, fill = saa)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.4) +
  geom_hline(yintercept = 7000) +
  geom_line(aes(y = mean, color = saa)) +
  scale_y_continuous("Protected koala habitat (ha)", limits = c(5000, 12000)) +
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  facet_wrap(~model) +
  ggpubr::theme_pubr()

ggsave("plots/saa_p1.png", saa_p1)
ggsave("plots/saa_p2.png", saa_p2)


obj_values_saa <- list()
obj_values_saa$no_recourse <- read_csv("data/saa_cost_nr.csv")
obj_values_saa$add_only <- read_csv("data/saa_cost_add.csv")
obj_values_saa$terminate_only <- read_csv("data/saa_cost_terminate.csv")
obj_values_saa$full_recourse <- read_csv("data/saa_cost_full.csv")
obj_values_saa_df <- obj_values_saa %>%
  lapply(function(x) {
    x %>% as.matrix()
    colnames(x) <- NA
    x <- as.matrix(x)
    c(x)
    }) %>%
  bind_rows(.id = 'model')
obj_values_saa_plot <- obj_values_saa_df %>%
  as.data.frame() %>%
  pivot_longer(colnames(obj_values_saa_df), names_to = 'model') %>%
  mutate(model = factor(model, colnames(obj_values_saa_df), c('None', 'Add', 'Terminate', 'Full'))) %>%
  ggplot(aes(x = value, y = model)) +
  geom_boxplot() +
  scale_x_continuous("Conservation cost", labels = scales::unit_format(prefix = "A$", suffix = "M",scale = 1e-6)) +
  scale_y_discrete('') +
  ggpubr::theme_pubr()
obj_values_saa_plot

ggsave("plots/obj_values_saa.png", obj_values_saa_plot)
