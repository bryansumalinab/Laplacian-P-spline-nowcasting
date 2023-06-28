Nowcasting.NB <- function(data, # Data containing date of occurence and reporting
                                date.now, # Nowcast date
                                max.delay, # Maximum delay
                                day.effect, # True if include day of the week effect
                                ref.day # reference category for day of the week effect
){
  
  Sys.setlocale("LC_TIME", "English") # set language to English
  
  # Data setup
  date.start <- as.Date(min(data$Date)) # Starting date
  T.now <- as.numeric(date.now - date.start) + 1
  
  data <- data %>% mutate_at(vars(t, d, Cases), as.numeric) %>% 
    mutate_at(vars(Date, Rep.date), as.Date) %>%
    filter(Date >= date.start & Date <= date.now) %>% 
    filter(d <= max.delay) %>%
    mutate(Reported = ifelse(t + d <= T.now, "Reported", "Not yet reported") 
           %>% factor(levels = c("Reported", 
                                 "Not yet reported")))
  
  # Hyperparameters for Gamma prior of delta
  a.delta <- b.delta <- 10^(-5)
  
  # Prior for overdispersion parameter phi
  a.phi <- b.phi <- 10^(-5)
  
  # nu prior parameter for the penalty
  nu <- 3
  
  D <- max(data$d) # Maximum delay
  TT <- max(data$t) # Maximum time
  
  # Rows for not yet reported
  nyr <- which(data$Reported == "Not yet reported")
  
  # Data for not yet reported (used for estimation)
  data_nyr <- data[nyr,]
  
  # Time and delay
  t <- unique(data$t)
  d <- unique(data$d)
  
  # Model matrices
  
  Kt = 40 # Number of B-splines for time dimension
  Kd = 10 # Number of B-splines for delay dimension
  
  # B-spline basis matrix
  Bt <- cubicbs(t,lower = min(t),upper = max(t),K = Kt)$Bmatrix
  Bd <- cubicbs(d,lower = min(d),upper = max(d),K = Kd)$Bmatrix
  
  # Two dimensional B-spline matrix
  B <- kronecker(Bd, Bt)
  
  y <- data$Cases[-nyr] # Reported cases
  
  # Difference order of the penalty
  penorder <- 2
  
  # Penalty for column (delay dimension)
  Dd <- diag(Kt)
  for (k in 1:penorder) Dd <- diff(Dd)
  Pd <- t(Dd) %*% Dd
  Pd <- Pd + diag(1e-12, Kt)
  
  # Penalty for row (time dimension)
  Dt <- diag(Kd)
  for (k in 1:penorder) Dt <- diff(Dt)
  Pt <- t(Dt) %*% Dt
  Pt <- Pt + diag(1e-12, Kd)
  
  
  if(day.effect == T){
    
    # Add day of the week
    data$Day <- weekdays(as.Date(data$Rep.date))
    data$Day <- relevel(factor(data$Day),ref.day)
    
    zeta <- 1e-05 # precision week effects coefficient
    
    # Design matrix for day of the week effect
    Z<-model.matrix(~ Day, data = data)
    
    # Global design matrix
    X1 <- cbind(B, Z)
    
    p <- ncol(Z)-1
    
    X_nyr <- X1[nyr,] # Design matrix for not yet reported cases 
    X <- X1[-nyr,] # # Design matrix for reported cases
    
    # Precision matrix for B-spline parameters
    Pv <- function(v) exp(v[1])*(kronecker(diag(1,Kd),Pd)) + 
      exp(v[2])*(kronecker(Pt,diag(1,Kt)))
    
    # Precision matrix for parameter xi
    Qv <- function(v) as.matrix(bdiag(diag(zeta, p+1), Pv(v)))
  } else{
    # Full design matrix
    X1 <- B
    
    X_nyr <- X1[nyr,] # Design matrix for not yet reported cases
    X <- X1[-nyr,] # Design matrix for reported cases
    
    # Precision matrix for B-spline parameters
    Pv <- function(v) exp(v[1])*(kronecker(diag(1,Kd),Pd)) + 
      exp(v[2])*(kronecker(Pt,diag(1,Kt)))
    
    # Precision matrix for parameter xi
    Qv <- Pv
  }
  
  # Negative binomial GLM with log-link
  mu.nb <- function(xi) exp(as.numeric(X %*% xi))
  var.nb <- function(xi, v) mu.nb(xi) + (1/exp(v[3]))*(mu.nb(xi)^2)
  W.nb <- function(xi, v) diag(((exp(as.numeric(X %*% xi)))^2)*(1/var.nb(xi, v)))
  D.nb <- function(xi) diag(1/mu.nb(xi))
  M.nb <- function(xi) diag(y - mu.nb(xi))
  V.nb <- function(xi, v) diag(mu.nb(xi) * (1/var.nb(xi, v) - 
                                              (mu.nb(xi)/(var.nb(xi, v)^2)) * 
                                              (1 + 2*mu.nb(xi)*(1/exp(v[3])))))
  gamma.nb <- function(xi, v) exp(v[3]) * log(mu.nb(xi)/(mu.nb(xi) + exp(v[3])))
  bgamma.nb <- function(xi, v) - (exp(v[3])^2) * log(exp(v[3])/(exp(v[3]) + mu.nb(xi)))
  
  # Log conditional posterior of xi given v
  log_pxi <- function(xi, v) {
    value <- (1/exp(v[3])) * sum((y * gamma.nb(xi, v)) - bgamma.nb(xi, v)) - .5 * t(xi) %*% Qv(v[1:2]) %*% xi
    return(value)
  }
  
  # Gradient of parameter xi
  Grad.logpxi <- function(xi,v){
    value <- t(X)%*%W.nb(xi,v)%*%D.nb(xi)%*%(y - mu.nb(xi)) - Qv(v[1:2])%*%xi
    as.numeric(value)
  }
  
  # Hessian of parameter xi
  Hess.logpxi <- function(xi,v){
    value <- t(X)%*%M.nb(xi)%*%V.nb(xi, v)%*%X - 
      t(X)%*%W.nb(xi, v)%*%X - Qv(v[1:2])
    value
  }
  
  # Initial values for log-penalty and log-overdispersion parameter
  v_init = c(1,1,1)
  
  # Initial estimate for xi
  xi_init <- NR_xi(xi0 = rep(0,dim(X)[2]), v = v_init, 
                   Hess.logpxi = Hess.logpxi, 
                   Grad.logpxi = Grad.logpxi, 
                   log_pxi = log_pxi)
  
  # Estimate of overdispersion parameter (phi) 
  n0 <- sum(y==0)
  phifun <- function(rho){
    val <- n0/length(y) - (1+mean(y)/rho)^(-rho)
    return(val)
  }
  vphi <- log(uniroot(f=phifun, lower = 1e-05, upper = 30)$root+2)
  
  # Log-conditional posterior of penalty on the time dimension
  log_pv1cond <- function(v){
    
    v <- c(v,0,vphi)
    
    e1 <- eigen(Pv(v[1:2]),only.values = T)$values
    e2 <- eigen(-t(X)%*%(M.nb(xi = xi_init)%*%V.nb(xi = xi_init, v = v)-W.nb(xi = xi_init, v = v))%*%X +
                  Qv(v[1:2]),only.values = T)$values
    
    value <- sum((1/exp(v[3]))*((y * gamma.nb(xi = xi_init, v)) - bgamma.nb(xi = xi_init, v)) +
                   lgamma(y + exp(v[3])) - lgamma(exp(v[3]))) +
      0.5*sum(sapply(e1[e1>0], log)) - 0.5 * sum((xi_init * Qv(v[1:2])) %*% xi_init) +
      a.phi*v[3] - b.phi*exp(v[3]) - 0.5*sum(sapply(e2[e2>0],log)) +
      0.5*nu*(v[1]+v[2]) - (0.5*nu + a.delta)*(log(b.delta + 0.5*nu*exp(v[1]))+log(b.delta + 0.5*nu*exp(v[2])))
    
    return(value)
  }
  
  v1dom <- seq(-0.5,3,length = 10)
  v1img <- sapply(v1dom, log_pv1cond)
  v1star <- v1dom[which(v1img==max(v1img))]
  
  # Conditional posterior mode of v
  v_mode <- c(v1star,0, vphi)
  
  # Mode a posteriori estimate for xi
  xi_mode <- NR_xi(xi0 = xi_init, v = v_mode,
                   Hess.logpxi = Hess.logpxi, 
                   Grad.logpxi = Grad.logpxi, 
                   log_pxi = log_pxi)
  
  ## Nowcast for not yet reported
  mu_nyr <- exp(X_nyr%*%xi_mode)
  nowcast <- data
  nowcast[nyr,"Cases"] <- mu_nyr
  
  ## Nowcasted incidence (reported + nowcast)
  nowcast.incidence <- nowcast %>% group_by(t) %>%
    summarize(y = ceiling(sum(Cases)))
  
  ## Prediction Interval
  
  pred.inter <- predinterval.NB(data = data,
                             date.start = date.start,
                             date.now = date.now,
                             X_nyr = X_nyr,
                             v_mode = v_mode,
                             xi_mode = xi_mode,
                             Qv = Qv)
  
  # Delay Density
  mu_hat <- exp(X1[,1:dim(B)[2]]%*%xi_mode[1:dim(B)[2]])
  nowcast2 <- data
  nowcast2[,"Cases"] <- mu_hat
  
  cases_matrix <- as.data.frame(matrix(nowcast2$Cases,nrow = TT, ncol = D+1, byrow = F))
  delaydist <- t(apply(cases_matrix, MARGIN = 1, function(i) i/sum(i)))
  data_delay <- data.frame("Date" = as.Date(data$Date),"Delay" = data$d,"density" = as.numeric(delaydist))
  
  output <- list(data = data,
                 nowcast.incidence = nowcast.incidence,
                 nowcast = pred.inter,
                 v_mode = v_mode,
                 xi_mode = xi_mode,
                 delay = data_delay,
                 date.now = date.now)
  output
}
