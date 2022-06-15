Nowcasting_sim<-function(data,date.now,Kt=40,Kd=10,v_init=c(1,1)){
  
  #prior for penalty 
  a=b= 10^(-5) 
  nu=3
  zeta=1e-05 #precision week effects coefficient
  
  D=max(data$d)
  T=max(data$t)
  
  #rows for not yet reported
  nyr=which(data$Reported=="Not yet reported")
  
  
  data_nyr=data[nyr,]
  data_rep=data[-nyr,]
  
  # Extract dimensions
  t <- unique(data$t)
  d <- unique(data$d)
  
  # Model matrices
  
  # B-spline basis matrix for smooth surface
  Bt<-cubicbs(t,lower = min(t),upper = max(t),K=Kt)$Bmatrix
  Bd<-cubicbs(d,lower = min(d),upper = max(d),K=Kd)$Bmatrix
  Bt<-Bt[,-Kt]
  Bd<-Bd[,-Kd]
  B <- kronecker(Bd, Bt)
  X1<-B
  
  Kt=dim(Bt)[2]
  Kd=dim(Bd)[2]
  
  X_nyr=X1[nyr,] 
  X=X1[-nyr,] #model matrix for reported cases
  y=data$Cases[-nyr] #reported cases
  
  # Penalty matrix
  penorder = 2
  #Penalty for column
  Dd <- diag(Kt)
  for (k in 1:penorder) Dd <- diff(Dd)
  Pd <- t(Dd) %*% Dd
  Pd <- Pd + diag(1e-12, Kt)
  
  #Penalty for row
  Dt <- diag(Kd)
  for (k in 1:penorder) Dt <- diff(Dt)
  Pt <- t(Dt) %*% Dt
  Pt <- Pt + diag(1e-12, Kd)
  
  #For poisson with log-link
  mu<- function(xi) exp(as.numeric(X %*% xi))
  W <- function(xi) diag(exp(as.numeric(X %*% xi)))
  s <- function(gam) exp(gam)
  s1 <- function(gam) exp(gam)
  s2 <- function(gam) exp(gam)
  
  Pv=function(v) exp(v[1])*(kronecker(diag(1,Kd),Pd)) + 
    exp(v[2])*(kronecker(Pt,diag(1,Kt)))
  Qv=Pv
  
  # Log conditional posterior of latent field given v
  logplat <- function(latent, v) {
    post <- sum(y * as.numeric(X %*% latent) -s(X %*% latent)) - .5 * t(latent) %*% Qv(v) %*% latent
    return(post)
  }
  
  #gradient of latent parameter xi
  grad.logplat=function(xi,v){
    value=t(X)%*%(y-mu(xi)) - Qv(v)%*%xi
    as.numeric(value)
  }
  
  #hessian of latent parameter xi
  Hess.logplat=function(xi,v){
    value= -t(X)%*%W(xi)%*%X - Qv(v)
    value
  }
  
  # Laplace approximation to conditional latent field
  Laplace <- function(latent0, v){
    
    epsilon <- 1e-05 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      dlat <- as.numeric((-1) * solve(Hess.logplat(latent0, v),
                                      grad.logplat(latent0, v)))
      lat.new <- latent0 + dlat
      step <- 1
      iter.halving <- 1
      logplat.current <- logplat(latent0, v)
      while (logplat(lat.new, v) <= logplat.current) {
        step <- step * .5
        lat.new <- latent0 + (step * dlat)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((lat.new - latent0) ^ 2))
      iter <- iter + 1
      latent0 <- lat.new
      if(dist < epsilon) break
    }
    
    latentstar <- latent0 
    return(latentstar)
  }
  
  v_init=rep(1,2)
  xi_hat=Laplace(rep(0,dim(X)[2]), v = v_init)
  
  XWX <- t(X)%*%W(xi_hat)%*%X
  Xxi <- X%*%xi_hat
  
  log_pv <- function(v){
    e1<-eigen(XWX + Qv(v),only.values = T)$values
    e2<-eigen(Pv(v),only.values = T)$values
    
    #a1 <- 0.5*log(det(XWX + Qv(v)))
    a1 <- 0.5*sum(sapply(e1,log))
    a2 <- sum(y*(Xxi))
    a3 <- sum(s(Xxi))
    a4 <- 0.5 * sum((xi_hat * Qv(v)) %*% xi_hat)
    #a5 <- 0.5*log(det(Pv(v)))
    a5 <- 0.5*sum(sapply(e2, log))
    a6 <- (0.5*nu + a)*(log(b + 0.5*nu*exp(v[1]))+log(b + 0.5*nu*exp(v[2])))
    a7 <- 0.5*nu*(v[1]+v[2])
    
    value=-a1+a2-a3-a4+a5-a6+a7
    return(value)
  }
  
  # Mode of v
  NR_v <- function(v0){
    
    epsilon <- 1e-05 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    for (k in 1:maxiter) {
      dlat <- as.numeric((-1) * solve(hessian(func=log_pv,x=v0),
                                      grad(func=log_pv,x=v0,method = "simple")))
      v.new <- v0 + dlat
      step <- 1
      iter.halving <- 1
      logpv.current <- log_pv(v0)
      while (log_pv(v.new) <= logpv.current) {
        step <- step * .5
        v.new <- v0 + (step * dlat)
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
  
  v_mode=NR_v(v0=v_init) 
  v_mode
  
  xi_mode=Laplace(xi_hat, v = v_mode)
  
  ##Nowcasting
  mu_nyr=exp(X_nyr%*%xi_mode) #nowcast
  nowcast<-data
  nowcast[nyr,"Cases"]<-mu_nyr
  
  #Time that has not yet reported cases is t=T-(D-1),...,T
  t.now=T-(D-1)
  nowcast<-nowcast %>% filter(t>=t.now) %>% select(Date,Cases)
  
  
  data1 <- data %>%
    mutate(Reported = Reported %>%
             factor(
               levels = c("Reported", "Not yet reported", "Nowcast"),
               labels = c("Reported", "Not yet reported", "Nowcast"))) %>% 
    group_by(Date, Reported) %>%
    summarize(Cases = sum(Cases)) %>% 
    mutate(Date=as.Date(Date)) %>% as.data.frame() %>%
    filter(Date>date.now-(2*D))
  
  data2<-nowcast %>% group_by(Date) %>%  summarize(Cases = sum(Cases))%>% 
    mutate(Date=as.Date(Date))
  
  
  ####################################################
  ###Confidence Interval
  W1 <- function(xi) diag(exp(as.numeric(X_nyr %*% xi)))
  Q=Qv(v_mode)
  sigma_xi_mode=solve(t(X_nyr)%*%W1(xi_mode)%*%X_nyr + Q) #covariance for xi
  
  
  #Compute mean and variance for log(mu_t.d)
  logmu=c()
  logmu.var=c()
  for(i in 1:dim(X_nyr)[1]){
    logmu[i]=X_nyr[i,]%*%xi_mode
    logmu.var[i]= X_nyr[i,]%*%sigma_xi_mode%*%X_nyr[i,]
  }
  
  r.pois<-list()
  B=1000
  
  for (i in 1:length(logmu)) {
    rn=rnorm(B, mean =logmu[i], sd=sqrt(logmu.var[i]))
    mu=exp(rn)
    r.pois[[i]]=rpois(B,lambda=mu)
  }
  
  data_nyr$nowcast<-exp(logmu)
  data_CI=data_nyr %>% 
    group_by(Date) %>% 
    summarise(nowcast=sum(nowcast)) %>% 
    mutate(Date=as.Date(Date))
  
  
  rownames(data_nyr)<-NULL
  days=seq(date.now-(D-1),date.now,by="day") # dates that has nowcast
  
  #sum of poisson samples for each t (days)
  for (i in 1:dim(data_CI)[1]) {
    r.pois.sum=Reduce("+",r.pois[as.numeric(rownames(data_nyr[data_nyr$Date == days[i],]))])
    data_CI[i,"lower_nyr"]=quantile(r.pois.sum,probs = 0.025)
    data_CI[i,"upper_nyr"]=quantile(r.pois.sum,probs = 0.975)
  }
  
  CI_nyr<-data.frame(data2,data_CI[,c("nowcast","lower_nyr","upper_nyr")])
  
  data_rep_cases=data1 %>% filter(Date %in% CI_nyr$Date) %>% 
    filter(Reported=="Reported") %>% data.frame() %>% 
    add_row(Date = max(data$Date), Reported = "Reported", Cases = 0)
  #add row 0 cases for reported on nowcast date,
  #to match number of rows for next code
  
  data_cases=data1 %>% filter(Date %in% CI_nyr$Date) %>%
    group_by(Date) %>% summarise(ObsCases=sum(Cases))
  
  CI_nyr$Rep_Cases=data_rep_cases$Cases
  CI_nyr$ObsCases=data_cases$ObsCases
  
  result=CI_nyr %>% 
    mutate(lower=lower_nyr+Rep_Cases,
           upper=upper_nyr+Rep_Cases) %>%
    mutate(x=ifelse(ObsCases>=lower & ObsCases<=upper,1,0)) %>%
    mutate(bias=Cases-ObsCases,
           mape=abs((Cases-ObsCases)/ObsCases),
           smape=abs(Cases-ObsCases)/(((abs(Cases)+abs(ObsCases)))/2)) %>% 
    select(-c(nowcast,lower_nyr,upper_nyr,Rep_Cases))
  
  
  output<-list(nowcast=nowcast,
               v_mode=v_mode,
               xi_mode=xi_mode,
               result=result)
  output
}