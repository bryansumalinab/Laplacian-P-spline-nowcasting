# Laplace approximation to conditional posterior of xi
# using Newton-Raphson algorithm
NR_xi <- function(xi0, v, Hess.logpxi, Grad.logpxi, log_pxi){
  
  epsilon <- 1e-05 # Stop criterion
  maxiter <- 100   # Maximum iterations
  iter <- 0        # Iteration counter
  
  for (k in 1:maxiter) {
    dxi <- as.numeric((-1) * solve(Hess.logpxi(xi0, v),
                                   Grad.logpxi(xi0, v)))
    xi.new <- xi0 + dxi
    step <- 1
    iter.halving <- 1
    logpxi.current <- log_pxi(xi0, v)
    while (log_pxi(xi.new, v) <= logpxi.current) {
      step <- step * .5
      xi.new <- xi0 + (step * dxi)
      iter.halving <- iter.halving + 1
      if (iter.halving > 30) {
        break
      }
    }
    dist <- sqrt(sum((xi.new - xi0) ^ 2))
    iter <- iter + 1
    xi0 <- xi.new
    if(dist < epsilon) break
  }
  
  xistar <- xi0 
  return(xistar)
}
