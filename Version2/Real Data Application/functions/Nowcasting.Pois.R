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
  
  # # Prior for overdispersion parameter phi
  # a.phi <- b.phi <- 10^(-5)
  
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
  
  # # Laplace approximation to conditional posterior of parameter xi
  # Laplace <- function(latent0, v){
  #   
  #   epsilon <- 1e-05 # Stop criterion
  #   maxiter <- 100   # Maximum iterations
  #   iter <- 0        # Iteration counter
  #   
  #   for (k in 1:maxiter) {
  #     dlat <- as.numeric((-1) * solve(Hess.logplat(latent0, v),
  #                                     grad.logplat(latent0, v)))
  #     lat.new <- latent0 + dlat
  #     step <- 1
  #     iter.halving <- 1
  #     logplat.current <- logplat(latent0, v)
  #     while (logplat(lat.new, v) <= logplat.current) {
  #       step <- step * .5
  #       lat.new <- latent0 + (step * dlat)
  #       iter.halving <- iter.halving + 1
  #       if (iter.halving > 30) {
  #         break
  #       }
  #     }
  #     dist <- sqrt(sum((lat.new - latent0) ^ 2))
  #     iter <- iter + 1
  #     latent0 <- lat.new
  #     if(dist < epsilon) break
  #   }
  #   
  #   latentstar <- latent0 
  #   return(latentstar)
  # }
  
  # # Initial estimate for xi
  # xi_init <- Laplace(rep(0,dim(X)[2]), v = v_init)
  
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
  
  # # Newton-Raphson to find posterior mode of v
  # NR_v  <-  function(v0){
  #   
  #   epsilon <- 1e-05 # Stop criterion
  #   maxiter <- 100   # Maximum iterations
  #   iter <- 0        # Iteration counter
  #   for (k in 1:maxiter) {
  #     dlat <- as.numeric((-1) * solve(hessian(func = log_pv,x = v0),
  #                                     grad(func = log_pv,x = v0,method = "simple")))
  #     v.new <- v0 + dlat
  #     step <- 1
  #     iter.halving <- 1
  #     logpv.current <- log_pv(v0)
  #     while (log_pv(v.new) <= logpv.current) {
  #       step <- step * .5
  #       v.new <- v0 + (step * dlat)
  #       iter.halving <- iter.halving + 1
  #       if (iter.halving > 30) {
  #         break
  #       }
  #     }
  #     dist <- sqrt(sum((v.new - v0) ^ 2))
  #     iter <- iter + 1
  #     v0 <- v.new
  #     if(dist < epsilon) break
  #   }
  #   
  #   vstar <- v0
  #   return(vstar)
  # }
  
  # # Mode of v
  # v_mode <- NR_v(v0 = v_init)
  # Posterior mode of v
  v_mode <- NR_v(v0 = v_init, log_pv = log_pv)
  
  # # Mode a posteriori estimate for xi
  # xi_mode <- Laplace(xi_init, v = v_mode)
  
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
  
  # # Time that has not yet reported cases is t = T-(D-1),...,T
  # t.now <- T-(D-1)
  # nowcast <- nowcast %>% filter(t >= t.now) %>% select(Date,Cases)
  # 
  # # Summarize the data and nowcast results 
  # # to be used for summarizing prediction interval in the next codes
  # data1 <- data %>%
  #   mutate(Reported = Reported %>%
  #            factor(
  #              levels = c("Reported", "Not yet reported", "Nowcast"),
  #              labels = c("Reported", "Not yet reported", "Nowcast"))) %>% 
  #   group_by(Date, Reported) %>%
  #   summarize(Cases = sum(Cases), .groups = "drop") %>% 
  #   mutate(Date = as.Date(Date)) %>% as.data.frame() %>%
  #   filter(Date>date.now-(2*D))
  # 
  # data2 <- nowcast %>% group_by(Date) %>%  summarize(Cases = sum(Cases))%>% 
  #   mutate(Date = as.Date(Date))
  
  ## Prediction Interval
  
  pred.inter <- predinterval.pois(data = data,
                             date.start = date.start,
                             date.now = date.now,
                             X_nyr = X_nyr,
                             v_mode = v_mode,
                             xi_mode = xi_mode,
                             Qv = Qv)
  
  # ####################################################
  # 
  # ### Confidence Interval
  # W1 <- function(xi) diag(exp(as.numeric(X_nyr %*% xi)))
  # Q <- Qv(v_mode)
  # sigma_xi_mode <- solve(t(X_nyr)%*%W1(xi_mode)%*%X_nyr + Q) # covariance for xi
  # 
  # 
  # # Compute mean and variance for log(mu_t.d)
  # logmu <- c()
  # logmu.var <- c()
  # for(i in 1:dim(X_nyr)[1]){
  #   logmu[i] <- X_nyr[i,]%*%xi_mode
  #   logmu.var[i] <- X_nyr[i,]%*%sigma_xi_mode%*%X_nyr[i,]
  # }
  # 
  # # Generate Poisson samples
  # r.pois <- list()
  # B <- 1000
  # for (i in 1:length(logmu)) {
  #   rn <- rnorm(B, mean = logmu[i], sd = sqrt(logmu.var[i]))
  #   mu <- exp(rn)
  #   r.pois[[i]] <- rpois(B,lambda = mu)
  # }
  # 
  # data_nyr$nowcast <- exp(logmu)
  # data_CI <- data_nyr %>% 
  #   group_by(Date) %>% 
  #   summarise(nowcast = sum(nowcast)) %>% 
  #   mutate(Date = as.Date(Date))
  # 
  # 
  # rownames(data_nyr) <- NULL
  # days <- seq(date.now-(D-1),date.now,by = "day") # dates that has nowcast
  # 
  # # Sum of poisson samples for each t (days)
  # for (i in 1:dim(data_CI)[1]) {
  #   r.pois.sum <- Reduce("+",r.pois[as.numeric(rownames(data_nyr[data_nyr$Date == days[i],]))])
  #   data_CI[i,"lower_nyr"] <- quantile(r.pois.sum,probs = 0.025)
  #   data_CI[i,"upper_nyr"] <- quantile(r.pois.sum,probs = 0.975)
  # }
  # 
  # CI_nyr <- data.frame(data2,data_CI[,c("nowcast","lower_nyr","upper_nyr")])
  # 
  # data_rep_cases <- data1 %>% filter(Date %in% CI_nyr$Date) %>% 
  #   filter(Reported == "Reported") %>% data.frame() %>% 
  #   add_row(Date = max(data$Date), Reported = "Reported", Cases = 0)
  # # Add row 0 cases for reported on nowcast date,
  # # to match number of rows for next code
  # 
  # data_cases <- data1 %>% filter(Date %in% CI_nyr$Date) %>%
  #   group_by(Date) %>% summarise(ObsCases = sum(Cases))
  # 
  # CI_nyr$Rep_Cases <- data_rep_cases$Cases
  # CI_nyr$ObsCases <- data_cases$ObsCases
  # 
  # result <- CI_nyr %>% 
  #   mutate(lower = lower_nyr+Rep_Cases,
  #          upper = upper_nyr+Rep_Cases) %>%
  #   mutate(x = ifelse(ObsCases >= lower & ObsCases <= upper,1,0)) %>%
  #   mutate(bias = Cases-ObsCases,
  #          mape = abs((Cases-ObsCases)/ObsCases),
  #          smape = abs(Cases-ObsCases)/(((abs(Cases)+abs(ObsCases)))/2),
  #          ciwdth = upper - lower) %>% 
  #   select(-c(nowcast,lower_nyr,upper_nyr,Rep_Cases))
  # 
  # 
  # output <- list(nowcast = nowcast,
  #                v_mode = v_mode,
  #                xi_mode = xi_mode,
  #                result = result)
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



