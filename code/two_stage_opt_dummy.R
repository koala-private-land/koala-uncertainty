source('code/fcn_two_stage_opt.R')

# Create dummy data ------------

# Number of properties
n <- 1000

# Area of property
h <- runif(n) * 10

# Cost of conservation in two time periods
c <- h + rnorm(n, 1, 1)
c <- matrix(rep(c,2), ncol = 2)

# Probability of each scenario (18 climate models)
s = 18
p = rep(1,s)/s

# Koala habitat (quality metric) under each climate realisation, in each timestep
t = 7
a_current <- h + rnorm(n, 0, 0.8)
a_matrix <- matrix(rep(a_current, t), ncol = t)
fcn_clim_data_sim <- function(a_matrix) {
  t_move <- stats::arima.sim(n = t, model = list(ar = 0.8))
  t_move <- t_move - t_move[1]
  a <- a_matrix + matrix(rep(t_move, n), ncol = t) + cbind(rep(0, n), matrix(rnorm(n*(t-1)), ncol = t-1))
}
a <- replicate(s, fcn_clim_data_sim(a_matrix), simplify = F)

# Recourse costs
beta <- 0.2 # Cost of adding a new property for conservation, per hectare
gamma <- 0.1 # Cost of terminating the contract of a property at year 2055

results <- fcn_two_stage_opt(cost = c, quality = a, area = h, target = c(80, 80, 80, 80, 100, 100, 100),
                             t_prime = 4, type = 'recourse')
results$x[1:n] %*% h
