library(gurobi)
fcn_lp_opt <- function(obj, A, b) {
  model <- list()
  model$A <- A
  model$obj <- obj
  model$rhs <- b
  model$sense <- '<'
  model$vtype <- 'C'
  model$lb <- 0
  model$ub <- 1
  model$modelsense <- 'min'
  params <- list()
  #params$NonConvex <- 2
  result <- gurobi(model, params)
  result
}

fcn_mv_opt <- function(obj, A, b, covmat, lambda) {
  model <- list()
  model$A <- A
  model$Q <- covmat * lambda 
  model$obj <- obj * (1-lambda)
  model$rhs <- b
  model$sense <- '<'
  model$vtype <- 'C'
  model$lb <- 0
  model$ub <- 1
  model$modelsense <- 'min'
  params <- list()
  params$NonConvex <- 2
  #params$NonConvex <- 2
  result <- gurobi(model, params)
}
