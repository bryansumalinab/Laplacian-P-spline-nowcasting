predinterval.pois <- function(data, # data
                         date.start, # Starting date
                         date.now, # Nowcast date
                         X_nyr, # Design matrix for not yet reported cases
                         v_mode, # Posterior mode of v
                         xi_mode, # Mode a posteriori estimate for xi
                         Qv # Precision matrix for parameter xi
){
  # Rows for not yet reported
  nyr <- which(data$Reported == "Not yet reported")
  # Data for not yet reported (used for estimation)
  data_nyr <- data[nyr,]
  y <- data_nyr$Cases # Reported cases
  
  ## Nowcasting
  mu_nyr <- exp(X_nyr%*%xi_mode) # nowcast
  nowcast <- data
  nowcast[nyr,"Cases"] <- mu_nyr
  
  # Time that has not yet reported cases is t = T-(D-1),...,T
  D <- max(data$d) # Maximum delay
  TT <- max(data$t) # Maximum time
  t.now <- max(data$t)-(D-1)
  nowcast.nyr <- nowcast %>% filter(t >= t.now) %>% select(Date,Cases)
  
  # Summarize the data and nowcast results
  # to be used for summarizing prediction interval in the next codes
  data1 <- data %>%
    mutate(Reported = Reported %>%
             factor(
               levels = c("Reported", "Not yet reported", "Nowcast"),
               labels = c("Reported", "Not yet reported", "Nowcast"))) %>%
    group_by(Date, Reported) %>%
    summarize(Cases = sum(Cases), .groups = "drop") %>%
    mutate(Date = as.Date(Date)) %>% as.data.frame() %>%
    filter(Date>date.now-(2*D))
  
  data2 <- nowcast.nyr %>% group_by(Date) %>%  summarize(Cases = sum(Cases))%>%
    mutate(Date = as.Date(Date))
  
  # Covariance of xi
  ### Confidence Interval
  W.po.nyr <- function(xi) diag(exp(as.numeric(X_nyr %*% xi)))
  Qv.nyr <- Qv(v_mode)
  sigma_xi <- solve(t(X_nyr)%*%W.po.nyr(xi_mode)%*%X_nyr + Qv.nyr) # covariance for xi
  
  # Compute mean and variance for log(mu_t.d)
  logmu <- c()
  logmu.var <- c()
  for(i in 1:dim(X_nyr)[1]){
    logmu[i] <- X_nyr[i,]%*%xi_mode
    logmu.var[i] <- X_nyr[i,]%*%sigma_xi%*%X_nyr[i,]
  }
  
  # # Generate negative binomial samples
  # r.nb <- list()
  # N <- 1000
  # for (i in 1:length(logmu)) {
  #   rn <- rnorm(N, mean = logmu[i], sd = sqrt(logmu.var[i]))
  #   mu <- exp(rn)
  #   r.nb[[i]] <- MASS::rnegbin(n = N, mu = mu, theta = exp(v_mode[3]))
  # }
  
  set.seed(12345)
  # Generate Poisson samples
  r.pois <- list()
  N <- 1000
  for (i in 1:length(logmu)) {
    rn <- rnorm(N, mean = logmu[i], sd = sqrt(logmu.var[i]))
    mu <- exp(rn)
    r.pois[[i]] <- rpois(N,lambda = mu)
  }
  
  data_nyr$nowcast <- exp(logmu)
  data_CI <- data_nyr %>% 
    group_by(t) %>% 
    summarise(nowcast = sum(nowcast))
  
  T.now <- as.numeric(date.now - date.start) + 1 
  days <- seq(T.now-(D-1),T.now) # dates that has nowcast
  
  # Sum of negative binomial samples for each t (days)
  for (i in 1:dim(data_CI)[1]) {
    r.pois.sum <- Reduce("+",r.pois[which(data_nyr$t == days[i])])
    data_CI[i,"lower_nyr"] <- quantile(r.pois.sum,probs = 0.025)
    data_CI[i,"upper_nyr"] <- quantile(r.pois.sum,probs = 0.975)
  }
  
  CI_nyr <- data.frame(data2,data_CI[,c("t","nowcast","lower_nyr","upper_nyr")])
  
  data_rep_cases <- data1 %>% filter(Date %in% CI_nyr$Date) %>%
    filter(Reported == "Reported") %>% data.frame()
  
  CI_nyr$Rep_Cases <- data_rep_cases$Cases
  
  value <- CI_nyr %>%
    mutate(`CI 95% lower` = ceiling(lower_nyr+Rep_Cases),
           `CI 95% upper` = ceiling(upper_nyr+Rep_Cases),
           Cases = ceiling(Cases)) %>%
    select(Date,Cases,`CI 95% lower`, `CI 95% upper`)
  
  return(value)
}