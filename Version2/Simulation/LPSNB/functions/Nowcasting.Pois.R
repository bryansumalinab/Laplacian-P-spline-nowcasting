Nowcasting.Pois <- function(data, # Data containing date of occurence and reporting
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
    # mutate(t = rep(1:length(unique(Date)), max.delay + 1)) %>%
    mutate(Reported = ifelse(t + d <= T.now, "Reported", "Not yet reported") 
           %>% factor(levels = c("Reported", 
                                 "Not yet reported")))
  
  # Hyperparameters for Gamma prior of delta
  a <- b <- 10^(-5)
  
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
  
  ###########################################
  if(day.effect == T){
    
    # Add day of the week
    data$Day <- weekdays(as.Date(data$Rep.date))
    data$Day <- relevel(factor(data$Day),ref.day)
    
    zeta <- 1e-05 #precision week effects coefficient
    
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
  ############################################
  
  # For Poisson GLM with log-link
  mu <- function(xi) exp(as.numeric(X %*% xi))
  W <- function(xi) diag(exp(as.numeric(X %*% xi)))
  s <- function(gam) exp(gam)
  s1 <- function(gam) exp(gam)
  s2 <- function(gam) exp(gam)
  
  # Log conditional posterior of parameter xi given v
  log_pxi <- function(xi, v) {
    post <- sum(y * as.numeric(X %*% xi) -s(X %*% xi)) - .5 * t(xi) %*% Qv(v) %*% xi
    return(post)
  }
  
  # Gradient of parameter xi
  Grad.logpxi <- function(xi,v){
    value <- t(X)%*%(y-mu(xi)) - Qv(v)%*%xi
    as.numeric(value)
  }
  
  # Hessian of parameter xi
  Hess.logpxi <- function(xi,v){
    value <- -t(X)%*%W(xi)%*%X - Qv(v)
    value
  }
  

  # Initial values for log-penalty and log-overdispersion parameter
  v_init = c(1,1)
  
  # Initial estimate for xi
  xi_init <- NR_xi(xi0 = rep(0,dim(X)[2]), v = v_init, 
                   Hess.logpxi = Hess.logpxi, 
                   Grad.logpxi = Grad.logpxi, 
                   log_pxi = log_pxi)
  
  # Log-posterior for the penalty- vector
  
  XWX <- t(X)%*%W(xi_init)%*%X
  Xxi <- X%*%xi_init
  
  log_pv <- function(v){
    e1 <- eigen(XWX + Qv(v),only.values = T)$values
    e2 <- eigen(Pv(v),only.values = T)$values
    
    a1 <- 0.5*sum(sapply(e1,log))
    a2 <- sum(y*(Xxi))
    a3 <- sum(s(Xxi))
    a4 <- 0.5 * sum((xi_init * Qv(v)) %*% xi_init)
    a5 <- 0.5*sum(sapply(e2, log))
    a6 <- (0.5*nu + a)*(log(b + 0.5*nu*exp(v[1]))+log(b + 0.5*nu*exp(v[2])))
    a7 <- 0.5*nu*(v[1]+v[2])
    
    value <- -a1+a2-a3-a4+a5-a6+a7
    return(value)
  }
  
  # Posterior mode of v
  v_mode <- NR_v(v0 = v_init, log_pv = log_pv)
  
  # Mode a posteriori estimate for xi
  xi_mode <- NR_xi(xi0 = xi_init, v = v_mode,
                   Hess.logpxi = Hess.logpxi, 
                   Grad.logpxi = Grad.logpxi, 
                   log_pxi = log_pxi)
  
  ## Nowcasting
  mu_nyr <- exp(X_nyr%*%xi_mode) # nowcast
  nowcast <- data
  nowcast[nyr,"Cases"] <- mu_nyr
  
  ## Nowcasted incidence (reported + nowcast)
  nowcast.incidence <- nowcast %>% group_by(t) %>%
    summarize(y = ceiling(sum(Cases)))
  
  
  ## Prediction Interval
  
  pred.inter <- predinterval.pois(data = data,
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



