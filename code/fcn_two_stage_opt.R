# Create and solve two-stage stochastic optimization problem ------
fcn_two_stage_opt <- function(cost, quality, area, target, t_prime, beta = 0, gamma = 0, type = 'recourse') {
  library(slam)
  library(gurobi)
  
  # Solve a two-stage conservation decision problem where in timestep t', all uncertainties are realised 
  # and the decision-maker can modify conservation plans
  #
  # Author: Frankie Cho
  # Date: 22 December 2022
  
  # Input arguments:
  # cost: a n x 2 matrix of conservation costs before and after the realisation of uncertainty
  # quality: a list with s elements, each element containing a n x t matrix of the quality metric for these properties for each timestep (e.g. suitability-weighted koala habitat)
  # area: vector of length n representing the area of each decision unit
  # target: a vector of length t 
  # t_prime (t'): Timestep when conservation plans can be modified and the climate scenario is known 
  # beta: per-area cost for adding new decision unit (on top of base cost) in the t' timestep (default = 0)
  # gamma: per-area cost for withdrawing a decision unit in the t' timestep (default = 0)
  # type: extent which recourse in conservation actions are allowed, with the following options:
  #   (1) "recourse": allows adding and removing decision units at timestep t' (default)
  #   (2) "only_new": only new decision units can be added at t', decision units selected before t' needs to continue to be selected
  #   (3) "only_remove": only 
  #
  # where n: number of decision units, s: number of scenarios, t: number of timesteps (derived within the function)
  #
  # Output:
  # A Gurobi object containing the solution of the problem
  c <- cost
  a <- quality
  h <- area
  K <- target
  n <- dim(c)[1]
  s <- length(a)
  t <- length(target)


  # Assemble elements
  cost_x <- rowSums(c)
  cost_y <- h * beta + c[,2]
  cost_w <- h * gamma - c[,2]
  
  mat_width <- n+n*s+n*s
  
  ## Constraint 1: habitat quality metric must reach K under all climate realisations before the timestep which the climate is known
  B1 <- rep(0, s*(t_prime-1)) # Negative because it is a > constraint
  ind <- 0
  A1 <- list()
  for (ss in 1:s) {
    for (tt in 1:(t_prime-1)) {
      ind <- ind + 1
      A1[[ind]] <- c(-t(a[[ss]][,tt]), rep(0, n*s+n*s))
      B1[ind] <- -K[tt]
    }
  }
  A1 <- do.call(rbind, A1)
  print('A1')
  
  ## Constraint 2: Koala habitat metric must reach K across all realisations after recourse actions
  ind <- 0
  B2 <- rep(0, s*(t-t_prime+1)) # Negative because it is a > constraint
  A2 <- list()
  for (ss in 1:s) {
    for (tt in t_prime:t) {
      ind <- ind + 1
      area_t_s <- t(a[[ss]][,tt])
      A2[[ind]] <- c(-area_t_s, -rep(area_t_s, s), rep(area_t_s, s))
      B2[ind] <- -K[tt]
    }
  }
  A2 <- do.call(rbind, A2)
  print('A2')
  
  ## Constraint 3: Additional properties (y) can only be activated if not selected in the first stage
  ind <- 0
  B3 <- rep(1, n*s)
  A3 <- list()
  for (ss in 1:s) {
    for (ii in 1:n) {
      ind <- ind + 1
      A3[[ind]] <- rep(0, mat_width)
      A3[[ind]][ii] <- 1
      A3[[ind]][n+ind] <- 1
    }
  }
  A3 <- do.call(rbind, A3)
  
  ## Vectorised fast version for constructing matrices -- WIP
  # nz <- 2 # nonzeros per row
  # ind <- 1:(s*n)
  # ii <- rep(1:n, s)
  # A3_i <- c(ii, n+ind)
  # A3_j <- c(ind, ind)
  # A3_v <- rep(1, s*n*nz)
  # A3 <- slam::simple_triplet_matrix(i=A3_i, j=A3_j, v=A3_v, ncol=mat_width, nrow = s*n)
  print('A3')
  
  ## Constraint 4: Properties can only be discontinued if covenant is signed in the first stage
  ind <- 0
  B4 <- rep(0, n*s)
  A4 <- list()
  for (ss in 1:s) {
    for (ii in 1:n) {
      ind <- ind + 1
      A4[[ind]] <- rep(0, mat_width)
      A4[[ind]][ii] <- -1
      A4[[ind]][n+n*s+ind] <- 1
    }
  }
  A4 <- do.call(rbind, A4)
  print('A4')
  
  # Solving the model
  model <- list()
  model$obj <- c(cost_x, rep(cost_y, s) * rep(p, each = n), rep(cost_w, s) * rep(p, each = n))
  model$modelsense <- 'min'
  model$A <- rbind(A1,A2,A3,A4)
  model$rhs <- c(B1,B2,B3,B4)
  model$sense <- '<'
  model$vtype <- 'C'
  model$lb <- rep(0, mat_width)
  model$ub <- rep(1, mat_width)
  
  params <- list()

  # Restricted model without recourse
  model_no_recourse <- model
  model_no_recourse$ub[(n+1):mat_width] <- 0
  
  # Only add new properties
  model_only_new <- model
  model_only_new$ub[(n+n*s+1):mat_width] <- 0
  
  # Only remove decision units
  model_only_remove <- model
  model_only_remove$ub[(n+1):(n+n*s)] <- 0
  
  # Select the appropriate model to solve
  model_solve <- switch(type, "recourse" = model, "only_new" = model_only_new, 
                        "only_remove" = model_only_remove, "no_recourse" = model_no_recourse)
  
  result <- gurobi(model_solve, params)
  return(result)
}

