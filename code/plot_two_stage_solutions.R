
properties_subset       <- readRDS('data/properties_subset.RData')

nsw_lga <- st_read( paste0(rdm_path, "\\planning_units\\planning_units.gdb"), layer = 'nsw_lga_pu_all') %>%
  st_transform(4326) %>%
  rename(lga = NSW_LGA__2) %>%
  st_simplify(preserveTopology = T, dTolerance = 2000)

solution_list <- list()
solution_list$no_recourse    <- read_csv('data/solution_no_recourse.csv')
solution_list$add_only       <- read_csv('data/solution_add_only.csv')
solution_list$terminate_only <- read_csv('data/solution_terminate_only.csv')
solution_list$full_recourse <- read_csv('data/solution_full_recourse.csv')
factor_recourse_var <- function(r) factor(r, c("no_recourse", "add_only", "terminate_only", "full_recourse"),
       c("No recourse", "Add only", "Terminate only", "Full recourse"))
solution <- solution_list %>% 
  bind_rows(.id = "recourse") %>%
  mutate(recourse = factor_recourse_var(recourse))

properties_solution <- properties_subset %>%
  right_join(solution)

obj_values <- read_csv("data/obj_values.csv")

obj_values_plot <- obj_values %>%
  mutate(Solution = factor(Solution, c("no_recourse", "add_only", "terminate_only", "full_recourse"),
                           c("No recourse", "Add only", "Terminate only", "Full recourse"))) %>%
  ggplot(aes(x = `Objective Value`, y = fct_rev(Solution))) +
  geom_bar(stat="identity", width = 0.5) +
  scale_y_discrete('') +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  ggpubr::theme_pubclean()
obj_values_plot

first_stage <- ggplot(properties_solution) +
  geom_sf(data = nsw_lga, fill = 'white') +
  geom_sf(data = properties_solution %>% filter(x == 0), color = 'gray60', size = 0.1) +
  geom_sf(data = properties_solution %>% filter(x == 1), color = '#3C5488FF') +
  facet_wrap(~recourse, nrow = 1) +
  theme_void()

properties_solution_long <- properties_solution %>%
  pivot_longer(cols = c("sum_y", "sum_w"), names_to = "variable", values_to = "decision") %>%
  mutate(variable = factor(variable, c("sum_y", "sum_w"), c("Add", "Terminate")))

second_stage <- ggplot(properties_solution_long) +
  geom_sf(data = nsw_lga, fill = 'white') +
  geom_sf(data = properties_solution_long %>% filter(decision < 1), color = 'gray60', size = 0.1) +
  geom_sf(data = properties_solution_long %>% filter(decision >= 1), aes(color = decision)) +
  scale_color_viridis_b() +
  facet_grid(rows = vars(variable), cols = vars(recourse)) +
  theme_void()

ggsave("plots/obj_values.png", obj_values_plot)
ggsave("plots/first_stage_solution.png", first_stage)
ggsave("plots/second_stage_solution.png", second_stage)
