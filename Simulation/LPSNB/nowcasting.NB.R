nowcasting.NB <- function(data){
  
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
  # Prior for overdispersion parameter phi
  a.phi <- 1e-05
  b.phi <- 1e-05
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
  
  # Precision matrix for parameter xi
  Pv <- function(v) exp(v[1])*(kronecker(Matrix::Diagonal(n = Kd),Pd)) +
    exp(v[2])*(kronecker(Pt,Matrix::Diagonal(n = Kt)))
  Qv <- Pv
  
  # Negative binomial GLM with log-link
  mu.nb <- function(xi) exp(as.numeric(X %*% xi))
  var.nb <- function(xi, v) {
    muval <- mu.nb(xi)
    res <- muval + (1 / exp(v[3])) * (muval ^ 2)
    return(res)
  }
  W.nb <- function(xi, v) {
    muval <- exp(as.numeric(X %*% xi))
    varval <- muval + (1 / exp(v[3])) * (muval ^ 2)
    res <- Matrix::Diagonal(x = ((muval) ^ 2) * (1 / varval))
    return(res)
  }
  
  D.nb <- function(xi) Matrix::Diagonal(x = 1/mu.nb(xi))
  M.nb <- function(xi) Matrix::Diagonal(x = y - mu.nb(xi))
  V.nb <- function(xi, v){
    muval <- mu.nb(xi)
    varval <-  muval + (1 / exp(v[3])) * (muval ^ 2)
    res <- Matrix::Diagonal(x = muval * (1/varval - (muval/(varval^2)) * 
                                           (1 + 2*muval*(1/exp(v[3])))))
    return(res)
  }
  gamma.nb <- function(xi, v) {
    muval <- mu.nb(xi)
    res <- exp(v[3]) * log(muval / (muval + exp(v[3])))
    return(res)
  }
  bgamma.nb <- function(xi, v) - (exp(v[3])^2) *
    log(exp(v[3])/(exp(v[3]) + mu.nb(xi)))
  
  # Log conditional posterior of xi given v
  log_pxi <- function(xi, v) {
    value <- (1/exp(v[3])) * sum((y * gamma.nb(xi, v)) - bgamma.nb(xi, v)) -
      .5 * t(xi) %*% Qv(v[1:2]) %*% xi
    return(as.numeric(value))
  }
  
  Grad.logpxi <- function(xi,v){
    muval <- exp(as.numeric(X %*% xi))
    varval <- muval + (1 / exp(v[3])) * (muval ^ 2)
    W.nbval <- Matrix::Diagonal(x = ((muval) ^ 2) * (1 / varval))
    D.nbval <- Matrix::Diagonal(x = 1/muval)
    value <- as.numeric(t(X)%*%W.nbval%*%D.nbval%*%(y - muval) -
                          Qv(v[1:2])%*%xi)
    return(value)
  }
  
  # Hessian of parameter xi
  Hess.logpxi <- function(xi,v){
    value <- t(X)%*%M.nb(xi)%*%V.nb(xi, v)%*%X -
      t(X)%*%W.nb(xi, v)%*%X - Qv(v[1:2])
    value
  }
  
  # Initial values for log-penalty and log-overdispersion parameter
  v_init <- c(1, 1, 1)
  
  NR_xi <- function(xi0, v){
    
    epsilon <- 1e-03 # Stop criterion
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
  
  xi_init <- NR_xi(xi0 = rep(0,dim(X)[2]), v = v_init)
  M.nbxi_init <- M.nb(xi = xi_init)
  muval_init <- mu.nb(xi_init)
  
  log_pv <- function(v){
    varval <-  muval_init + (1 / exp(v[3])) * (muval_init ^ 2)
    Vnb <- Matrix::Diagonal(x = muval_init * (1/varval - (muval_init/(varval^2)) *
                                                (1 + 2*muval_init*(1/exp(v[3])))))
    Wnb <- Matrix::Diagonal(x = ((muval_init) ^ 2) * (1 / varval))
    gammanb <- exp(v[3]) * log(muval_init / (muval_init + exp(v[3])))
    bgammanb <- (-1) * (exp(v[3])^2) * log(exp(v[3])/(exp(v[3]) + muval_init))
    
    e1 <- eigen(Pv(v[1:2]),only.values = T)$values
    e2 <- eigen(-t(X)%*%(M.nbxi_init%*%Vnb-Wnb)%*%X + Qv(v[1:2]),
                only.values = T)$values
    
    value <- sum((1/exp(v[3]))*((y * gammanb) - bgammanb) +
                   lgamma(y + exp(v[3])) - lgamma(exp(v[3]))) +
      0.5*sum(sapply(e1[e1>0], log)) -
      0.5 * sum((xi_init * Qv(v[1:2])) %*% xi_init) +
      a.phi*v[3] - b.phi*exp(v[3]) - 0.5*sum(sapply(e2[e2>0],log)) +
      0.5*nu*(v[1]+v[2]) - (0.5*nu + a.delta)*(log(b.delta + 0.5*nu*exp(v[1]))+
                                                 log(b.delta + 0.5*nu*exp(v[2])))
    
    return(as.numeric(value))
  }
  
  v_mode <- stats::optim(par = c(1,1,1), fn = log_pv,method = "Nelder-Mead",
                        control = list(fnscale = -1, reltol = 1e-5))$par
  
  # Conditional posterior mode of v
  xi_mode <- NR_xi(xi0 = xi_init, v = v_mode)
  
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
  sigma_xi <- -solve(Hess.logpxi(xi = xi_mode, v = v_mode))
  
  logmu <- numeric(dim(X_nyr)[1])
  logmu.var <- numeric(dim(X_nyr)[1])
  for(i in 1:dim(X_nyr)[1]){
    logmu[i] <- X_nyr[i,]%*%xi_mode
    logmu.var[i] <- X_nyr[i,]%*%sigma_xi%*%X_nyr[i,]
  }
  
  # Generate negative binomial samples
  r.nb <- list()
  N <- 1000
  for (i in 1:length(logmu)) {
    rn <- stats::rnorm(N, mean = logmu[i], sd = sqrt(logmu.var[i]))
    mu <- exp(rn)
    r.nb[[i]] <- stats::rnbinom(n = N, mu = mu, size = exp(v_mode[3]))
  }
  
  # Data for not yet reported cases
  data_nyr <- data[nyr,]
  data_nyr$nowcast <- exp(logmu)
  data_CI <- stats::aggregate(nowcast ~ t, data = data_nyr, FUN = sum)
  
  T.now <- as.numeric(date.now - date.start) + 1
  days <- seq(T.now-(D-1),T.now) # dates that have nowcast
  
  # Sum of negative binomial samples for each t (days)
  for (i in 1:dim(data_CI)[1]) {
    r.nb.sum <- Reduce("+",r.nb[which(data_nyr$t == days[i])])
    data_CI[i,"lower_nyr"] <- stats::quantile(r.nb.sum,probs = 0.025)
    data_CI[i,"upper_nyr"] <- stats::quantile(r.nb.sum,probs = 0.975)
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
  
  phi_estim <- data.frame("phi" = exp(v_mode[3]),
                          "description" = c("overdispersion parameter"))
  
  outputlist <- list(data = data,
                     cases.now = cases.now,
                     delay = data_delay,
                     lambda_estim = lambda_estim,
                     phi_estim = phi_estim)
  
  outputlist
}