# VOI under risk aversion
library(tidyverse)
library(mvtnorm)

# Value function across all states

## Baseline setup -- 2 actions and two states
n_states <- 2
n_actions <- 2
a1 <- c(1, 1)
a2 <- c(2, 0.9)
p <- c(0.5, 0.5)
action_state <- matrix(c(a1, a2), byrow = T, ncol = length(a1))

## Several states
#set.seed(1111)
#n_states <- 50
#n_actions <- 50
#A <- lhs::randomLHS(n = n_actions, k = n_states)
#action_state <- A
#p <- runif(n_states)
#p <- p/sum(p)

colnames(action_state) <- paste0('s', 1:n_states)
rownames(action_state) <- paste0('a', 1:n_actions)

gamma_seq <- seq(-30,30,.1)


# Utility function
CRRA <- function(c, gamma = 1) {
  if (gamma == 1) {
    return(log(c))
  } else if (gamma > 0) {
    return(c^(1-gamma) / (1-gamma))
  } else {
    return(NA)
  }
}

CRRA_inv <- function(c, gamma = 1) {
  if (gamma == 1) {
    return(exp(c))
  } else if (gamma > 0) {
    return((1-gamma)^(1/(1-gamma))*(c^(1/(1-gamma))))
  } else {
    return(NA)
  }
}

CARA <- function(c, alpha = 0) {
  if (alpha == 0 | is.na(alpha)) {
    return(c)
  } else {
    return((1-exp(-alpha*c))/alpha)
  }
}

CARA_inv <- function(c, alpha = 0) {
  if (alpha == 0 | is.na(alpha)) {
    return(c)
  } else {
    return(-(log(1-alpha*c))/alpha)
  }
}

# Expected Utility of decision-making with information
fcn_EVPI <- function(gamma) {
  U <- function(c) CARA(c, gamma)
  U_inv <- function(c) CARA_inv(c, gamma)
  utility_table <- U(action_state)
  max_a <- apply(utility_table, 2, which.max)  # Actions that maximise utility in each state (indices)
  max_EU <- which.max(utility_table %*% p) # Action that maximises expected utility (index)
  VPI_value <- action_state[cbind(max_a, seq_along(max_a))] %*% p - action_state[max_EU,] %*% p
  VPI <- U_inv(utility_table[cbind(max_a, seq_along(max_a))] %*% p) - U_inv(utility_table[max_EU,] %*% p)
  list(VPI = VPI, VPI_value = VPI_value, max_EU = max_EU)
}

evpi_seq <- lapply(gamma_seq, fcn_EVPI) %>%
  lapply(function(x) c(VPI = x$VPI, VPI_value = x$VPI_value, max_EU = x$max_EU)) %>%
  bind_rows(.id = 'gamma')
evpi_seq$gamma <- gamma_seq
evpi_seq %>%
  ggplot(aes(y = VPI, x = gamma, color = factor(max_EU))) +
  geom_point() +
  theme_minimal()

