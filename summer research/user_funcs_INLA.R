## Define CAR Model
inla.rgeneric.CAR.model <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {
  
  #Internal function
  interpret.theta <- function() {
    return(
      list(prec = exp(theta[1L]),
           rho = 1 / (1 + exp(-theta[2L])))
    )
  }
  
  graph <- function(){
    require(Matrix)
    
    return(Diagonal(nrow(W), x = 1) + W)
  }
  
  Q <- function() {
    require(Matrix)
    
    param <- interpret.theta()
    
    return(param$prec * (Diagonal(nrow(W), x = rowSums(W)) - param$rho * W) )
  }
  
  mu <- function()
  {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    return(numeric(0))
    
  }
  
  log.prior <- function() {
    param = interpret.theta()
    
    res <- dgamma(param$prec, 1, 5e-05, log = TRUE) + log(param$prec) +
      log(1) + log(param$rho) + log(1 - param$rho) 
    
    return(res)
  }
  
  initial <- function() {
    return(c(0, 0))
  }
  
  quit <- function() {
    return(invisible())
  }
  
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}
CAR.model <- inla.rgeneric.define(inla.rgeneric.CAR.model, W = W)
