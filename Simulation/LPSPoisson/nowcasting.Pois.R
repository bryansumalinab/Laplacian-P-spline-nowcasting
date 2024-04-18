nowcasting.Pois <- function(data){
  
  if (!require("Matrix", character.only = TRUE)) {
    message(paste("Package Matrix", "is required but not installed."))
  }
  
  date.start <- as.Date(min(data$Date)) # Starting date
  date.now <- as.Date(max(data$Date))   # Nowcast date
  max.delay <- max(data$d)              # Maximum delay
  data <- data[order(data$d), ]
  
  # Hyperparameters for Gamma prior of delta
  a.delta <- 1e-05
  b.delta <- 1e-05
  nu <- 3             # prior parameter for the penalty
  D <- max(data$d)    # Maximum delay
  TT <- max(data$t)   # Maximum time
  nyr <- which(data$Reported == "Not yet reported") # Rows for not yet reported
  
  # Time and delay
  t <- unique(data$t)
  d <- unique(data$d)
  
  # Model matrices
  Kt <- 40 # Number of B-splines for time dimension
  Kd <- 5 # Number of B-splines for delay dimension
  
  # B-spline basis matrix
  Bt <- cubicbs(t,lower = min(t),upper = max(t),K = Kt)$Bmatrix
  Bd <- cubicbs(d,lower = min(d),upper = max(d),K = Kd)$Bmatrix
  B <- kronecker(Bd, Bt)  # Two dimensional B-spline matrix
  y <- data$Cases[-nyr]   # Reported cases
  penorder <- 2           # Penalty order
  
  # Penalty for column (delay dimension)
  Dd <- Matrix::Diagonal(n = Kt)
  for (k in 1:penorder) Dd <- Matrix::diff(Dd)
  Pd <- Matrix::t(Dd) %*% Dd
  Pd <- Pd + Matrix::Diagonal(n = Kt, x = 1e-12)
  
  # Penalty for row (time dimension)
  Dt <- Matrix::Diagonal(n = Kd)
  for (k in 1:penorder) Dt <- Matrix::diff(Dt)
  Pt <- Matrix::t(Dt) %*% Dt
  Pt <- Pt + Matrix::Diagonal(n = Kd, x = 1e-12)
  
  # Full design matrix
  X1 <- B
  X_nyr <- X1[nyr,] # Design matrix for not yet reported cases
  X <- X1[-nyr,]    # Design matrix for reported cases
  Pv <- function(v) exp(v[1])*(kronecker(Matrix::Diagonal(n = Kd),Pd)) +
    exp(v[2])*(kronecker(Pt,Matrix::Diagonal(n = Kt)))
  Qv <- Pv # Precision matrix for parameter xi
  
  ############################################
  
  # For Poisson GLM with log-link
  mu <- function(xi) exp(as.numeric(X %*% xi))
  W <- function(xi) Matrix::Diagonal(x = exp(as.numeric(X %*% xi)))
  s <- function(gam) exp(gam)
  
  log_pxi <- function(xi, Qv) {
    value <- sum(y * as.numeric(X %*% xi) -s(X %*% xi)) - .5 * t(xi) %*% Qv %*% xi
    return(as.numeric(value))
  }
  
  Grad.logpxi <- function(xi, Qv){
    value <- t(X)%*%(y-mu(xi)) - Qv%*%xi
    as.numeric(value)
  }
  
  Hess.logpxi <- function(xi,Qv){
    value <- -t(X)%*%W(xi)%*%X - Qv 
    value
  }
  
  # Laplace approximation to conditional posterior of xi
  # using Newton-Raphson algorithm
  NR_xi <- function(xi0, Qv){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      dxi <- as.numeric((-1) * solve(Hess.logpxi(xi0, Qv),
                                     Grad.logpxi(xi0, Qv)))
      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- log_pxi(xi0, Qv)
      while (log_pxi(xi.new, Qv) <= logpxi.current) {
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
  
  
  # Initial values for log-penalty and log-overdispersion parameter
  v_init = c(1,1)
  
  # Initial estimate for xi
  Qv_init <- Qv(v_init)
  xi_init <- NR_xi(xi0 = rep(0,dim(X)[2]), Qv = Qv_init)
  
  # Log-posterior for the penalty- vector
  XWX <- t(X)%*%W(xi_init)%*%X
  Xxi <- X%*%xi_init
  
  log_pv <- function(v){
    Qv <- Qv(v)
    e1 <- eigen(XWX + Qv,only.values = T)$values
    e2 <- eigen(Pv(v),only.values = T)$values
    
    a1 <- 0.5*sum(sapply(e1,log))
    a2 <- sum(y*(Xxi))
    a3 <- sum(s(Xxi))
    a4 <- 0.5 * sum((xi_init * Qv) %*% xi_init)
    a5 <- 0.5*sum(sapply(e2, log))
    a6 <- (0.5*nu + a.delta)*(log(b.delta + 0.5*nu*exp(v[1]))+log(b.delta + 0.5*nu*exp(v[2])))
    a7 <- 0.5*nu*(v[1]+v[2])
    
    value <- -a1+a2-a3-a4+a5-a6+a7
    return(value)
  }
  
  # Posterior mode of v
  v_mode <- optim(par = v_init, 
                  fn = log_pv, 
                  method="Nelder-Mead", 
                  control = list(fnscale = -1, reltol = 1e-05))$par
  
  # Mode a posteriori estimate for xi
  Qv_mode <- Qv(v_mode)
  xi_mode <- NR_xi(xi0 = xi_init, Qv = Qv_mode)
  
  # Nowcast for not yet reported
  mu_nyr <- exp(X_nyr%*%xi_mode)
  nowcast <- data
  nowcast[nyr,"Cases"] <- mu_nyr
  
  # Nowcasted cases (reported + nowcast)
  cases.now <- stats::aggregate(Cases ~ t, data = nowcast,
                                FUN = function(x) ceiling(sum(x)))
  colnames(cases.now) <- c("t", "y")
  cases.now$Date <- unique(nowcast$Date)
  cases.now$CI95L <- NA
  cases.now$CI95R <- NA
  
  ##### Prediction Interval
  # Time that has not yet reported cases is t = T-(D-1),...,T.
  Date <- data$Date
  Cases <- data$Cases
  Reported <- data$Reported
  
  t.now <- max(data$t)-(D-1)
  nowcast.nyr <- subset(nowcast, t >= t.now, select = c(Date, Cases))
  
  # Summarize the data and nowcast results to be used for summarizing prediction
  # interval in next lines of code
  data1 <- stats::aggregate(Cases ~ Date + Reported, data = data, FUN = sum)
  data1$Reported <- factor(data1$Reported,
                           levels = c("Reported", "Not yet reported",
                                      "Nowcast"),
                           labels = c("Reported", "Not yet reported",
                                      "Nowcast"))
  data1$Date <- as.Date(data1$Date)
  data1 <- as.data.frame(data1)
  
  data2 <- stats::aggregate(Cases ~ Date, data = nowcast.nyr, FUN = sum)
  data2$Date <- as.Date(data2$Date)
  data2 <- as.data.frame(data2)
  
  # Covariance of xi
  sigma_xi <- -solve(Hess.logpxi(xi = xi_mode, Qv = Qv_mode))
  
  logmu <- numeric(dim(X_nyr)[1])
  logmu.var <- numeric(dim(X_nyr)[1])
  for(i in 1:dim(X_nyr)[1]){
    logmu[i] <- X_nyr[i,]%*%xi_mode
    logmu.var[i] <- X_nyr[i,]%*%sigma_xi%*%X_nyr[i,]
  }
  
  # Generate negative binomial samples
  r.Poiss <- list()
  N <- 1000
  for (i in 1:length(logmu)) {
    rn <- stats::rnorm(N, mean = logmu[i], sd = sqrt(logmu.var[i]))
    mu <- exp(rn)
    r.Poiss[[i]] <- stats::rpois(n = N, lambda = mu)
  }
  
  # Data for not yet reported cases
  data_nyr <- data[nyr,]
  data_nyr$nowcast <- exp(logmu)
  data_CI <- stats::aggregate(nowcast ~ t, data = data_nyr, FUN = sum)
  
  T.now <- as.numeric(date.now - date.start) + 1
  days <- seq(T.now-(D-1),T.now) # dates that have nowcast
  
  # Sum of negative binomial samples for each t (days)
  for (i in 1:dim(data_CI)[1]) {
    r.Poiss.sum <- Reduce("+",r.Poiss[which(data_nyr$t == days[i])])
    data_CI[i,"lower_nyr"] <- stats::quantile(r.Poiss.sum,probs = 0.025)
    data_CI[i,"upper_nyr"] <- stats::quantile(r.Poiss.sum,probs = 0.975)
  }
  
  CI_nyr <- data.frame(data2,data_CI[,c("t","nowcast","lower_nyr","upper_nyr")])
  data_rep_cases <- subset(data1, Date %in% CI_nyr$Date & Reported == "Reported")
  data_rep_cases <- as.data.frame(data_rep_cases)
  CI_nyr$Rep_Cases <- data_rep_cases$Cases
  
  CI_nyr$CI95L <- ceiling(CI_nyr$lower_nyr + CI_nyr$Rep_Cases)
  CI_nyr$CI95R <- ceiling(CI_nyr$upper_nyr + CI_nyr$Rep_Cases)
  
  cases.now[(nrow(cases.now) - max.delay + 1):nrow(cases.now),
            c("CI95L", "CI95R")] <- CI_nyr[, c("CI95L", "CI95R")]
  cases.now$status <- ifelse(is.na(cases.now$CI95L), "observed", "nowcasted")
  
  # Delay density
  mu_hat <- exp(X1[,1:dim(B)[2]]%*%xi_mode[1:dim(B)[2]])
  nowcast2 <- data
  nowcast2[,"Cases"] <- mu_hat
  
  cases_matrix <- as.data.frame(matrix(nowcast2$Cases,nrow = TT, ncol = D+1,
                                       byrow = F))
  delaydist <- t(apply(cases_matrix, MARGIN = 1, function(i) i/sum(i)))
  data_delay <- data.frame("Date" = as.Date(data$Date),
                           "Delay" = data$d,"density" = as.numeric(delaydist))
  
  lambda_estim <- data.frame("lambda" = exp(v_mode[1:2]),
                             "description" = c(" lambda_t (penalty for time)",
                                               "lambda_d (penalty for delay)"))
  outputlist <- list(data = data,
                     cases.now = cases.now,
                     delay = data_delay,
                     lambda_estim = lambda_estim)
  outputlist
}



