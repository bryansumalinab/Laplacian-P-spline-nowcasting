library(blapsr)
library(Matrix)
library(numDeriv)
library(tidyverse)
library(RColorBrewer)

source("Nowcasting_sim.R")

#function f1
f <- function(a1,a2,x) exp(a1 + a2*sin((2*pi*x)/150))
mu.sim=f(a1=2,a2=1,x=1:365) #f11
# mu.sim=f(a1=3,a2=1,x=1:365) #f12

# #function f2
# f <- function(a1,a2,x) exp(a1 + a2*sin((2*pi*x)/150) +0.2*sqrt(x))
# mu.sim=f(a1=0,a2=0.4,x=1:365)

p=c(0.1,0.4,0.2,0.1,0.1,0.05,0.05)
N=1000 #number of simulations
max.delay=7 #maximum delay
date.now=as.Date('2021-03-31') #set nowcast date

CI=matrix(nrow=max.delay,ncol=N)
bias=matrix(nrow=max.delay,ncol=N)
mape=matrix(nrow=max.delay,ncol=N)
smape=matrix(nrow=max.delay,ncol=N)

set.seed(12345)
for (j in 1:N) {
  tryCatch({
    
    simcases=matrix(nrow=365,ncol = 7)
    for (i in 1:365) {
      simcases[i,]=rmultinom(n=1,size = rpois(n=1,mu.sim[i]), prob = p)
    }
    
    dates=seq(from=as.Date("2021/1/1"),to=as.Date("2021/12/31") ,by = "day")
    
    data=data.frame(t=rep(1:365,times=7),d=rep(1:7,each=365),
                    Cases=as.vector(simcases),Date=rep(dates,times=7))
    data = data %>% mutate(Rep.date=Date+d)
    date.start=as.Date(min(data$Date))
    T.now=as.numeric(date.now - date.start) + 1
    data=data %>% mutate(Reported = ifelse(t + d <= T.now, "Reported",
                                           ifelse(t > T.now,"Future",
                                                  "Not yet reported")) %>%
                           factor(levels = c("Reported", "Not yet reported", "Future")))
    data=data %>%
      filter(Reported != "Future") %>%
      filter(d<=max.delay) %>%
      filter(Date<=date.now) %>%
      mutate(Date=as.Date(Date),
             Rep.date=as.Date(Rep.date),
             Reported=Reported %>%
               factor(levels = c("Reported", "Not yet reported")
                      ,labels = c("Reported", "Not yet reported"))) %>%
      data.frame()
    
    my_nowcast<-Nowcasting_sim(data=data,date.now=date.now)
    
    CI[,j]<-my_nowcast$result$x
    bias[,j]<-my_nowcast$result$bias
    mape[,j]<-my_nowcast$result$mape
    smape[,j]<-my_nowcast$result$smape
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  print(j)
}



######## (a) CI COVERAGE #######
round(apply(CI, MARGIN = 1, mean,na.rm=T)*100,2)

######## (b) MSE #######
round(apply(bias^2, MARGIN = 1, mean,na.rm=T),2)

######## (c) MAPE #######
#replace infinity be NA
mape[sapply(mape, is.infinite)] <- NA
round(apply(mape, MARGIN = 1, mean,na.rm=T)*100,2)

######## (d) SMAPE #######
round(apply(smape, MARGIN = 1, mean,na.rm=T)*100,2)
