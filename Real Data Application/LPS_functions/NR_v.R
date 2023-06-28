# Newton-Raphson to find posterior mode of v
NR_v <- function(v0, log_pv){
  
  epsilon <- 1e-05 # Stop criterion
  maxiter <- 100   # Maximum iterations
  iter <- 0        # Iteration counter
  for (k in 1:maxiter) {
    dv <- as.numeric((-1) * solve(hessian(func = log_pv,x = v0),
                                  grad(func = log_pv,x = v0,method  =  "simple")))
    v.new <- v0 + dv
    step <- 1
    iter.halving <- 1
    logpv.current <- log_pv(v0)
    while (log_pv(v.new) <= logpv.current) {
      step <- step * .5
      v.new <- v0 + (step * dv)
      iter.halving <- iter.halving + 1
      if (iter.halving > 30) {
        break
      }
    }
    dist <- sqrt(sum((v.new - v0) ^ 2))
    iter <- iter + 1
    v0 <- v.new
    if(dist < epsilon) break
  }
  
  vstar <- v0
  return(vstar)
}
