library(ggplot2)
library(LaplacesDemon)
library(mvtnorm)

source("~/GitHub/uncertainty/code/voi_simulations.R")

gamma_seq <- seq(-30,30,.1)
## Example 1: 2 actions and two states
e1 <- list()
n_states <- 2
n_actions <- 2
a1 <- c(1, 1)
a2 <- c(2, 0.9)
e1$p <- c(0.5, 0.5)
action_state <- matrix(c(a1, a2), byrow = T, ncol = length(a1))
colnames(action_state) <- paste0('s', 1:n_states)
rownames(action_state) <- paste0('a', 1:n_actions)
e1$action_state <- action_state
fcn_VOI_simulation(e1, gamma_seq) %>%
  ggplot(aes(x = gamma, y = VPI, color = factor(max_EU))) +
  geom_point() +
  theme_minimal()

## Several states
set.seed(1111)
n_states <- 50
n_actions <- 50
A <- lhs::randomLHS(n = n_actions, k = n_states)
action_state <- A
p <- runif(n_states)
p <- rep(1, n_states)
p <- p/sum(p)
n_y <- 5 # Partial information experiment that resolves uncertainty to n_y groups of possible outcomes
Y <- sort(rep(1:n_y, 1+n_states/n_y)[1:n_states])

## Draw from an inverse Wishart distribution
set.seed(11111)
n_states <- 100
n_actions <- 20

nsims <- 50
mu <- runif(n_actions)
wish_sim <- lapply(1:nsims, function(i) {
  set.seed(i)
  Sigma <- LaplacesDemon::rinvwishart(n_actions+1, diag(n_actions)) / 10
  prob <- list()
  action_state <- t(rmvnorm(n_states, mean = mu, sigma = Sigma))
  colnames(action_state) <- paste0('s', 1:n_states)
  rownames(action_state) <- paste0('a', 1:n_actions)
  prob$action_state <- action_state
  p <- runif(n_states)
  prob$p <- p/sum(p)
  n_y <- 5
  prob$Y <- sort(rep(1:n_y, 1+n_states/n_y)[1:n_states])
  prob
})

wish_sim_results <- lapply(wish_sim, fcn_VOI_simulation, gamma_seq = gamma_seq)
names(wish_sim_results) <- 1:nsims
wish_sim_table <- wish_sim_results %>%
  #lapply(function(x) mutate(x, VPI = VPI - x$VPI[x$gamma == 0])) %>%
  lapply(function(x) x$VPI) %>%
  bind_rows()
wish_sim_table$gamma_seq <- gamma_seq
wish_sim_table %>%
  pivot_longer(1:nsims) %>%
  ggplot(aes(x = gamma_seq, y = value, group = name)) +
  geom_line() +
  theme_minimal()

