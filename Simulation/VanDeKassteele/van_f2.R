rm(list = ls())
# Load packages
library(Matrix)
library(tidyverse)

# Source functions
list.files(path = "van_functions", full.names = TRUE) %>% 
  walk(source)


# Function f2
a1 <- 3
a2 <- 2
f <- function(a1,a2,x) 50 + exp(a1 + a2*sin((2*pi*x)/150))
mu.sim <- f(a1  = a1, a2 = a2, x = 1:365)

# Different nowcast dates (every end of the month) from March to November
dates.nowcast <- as.Date(c('2021-03-31','2021-04-30','2021-05-31','2021-06-30',
                           '2021-07-31','2021-08-31','2021-09-30','2021-10-31','2021-11-30'))

start.date <- as.Date("2021-01-01") # Set start date
nowcast.date <- dates.nowcast[1] # Set nowcast date

p <- c(0.1,0.4,0.2,0.1,0.1,0.05,0.05) # Delay probabilities
N <- 500 # Number of simulations
max.delay <- 7 # Maximum delay
overdisp <- 10


PI_cov <- matrix(nrow = max.delay,ncol = N)
PI_width <- matrix(nrow = max.delay,ncol = N)
Bias_mu <- matrix(nrow = max.delay,ncol = N)
Relbias_mu <- matrix(nrow = max.delay,ncol = N)

set.seed(12345)
for (j in 1:N) {
  tryCatch({
    simcases <- matrix(nrow = 365,ncol  =  7)
    for (i in 1:365) {
      # Generate multinomial samples as explained in the simulation procedure
      simcases[i,] <- rmultinom(n = 1,size = MASS::rnegbin(n = 1, mu = mu.sim[i], theta = overdisp), prob = p)
    }
    
    # Create a vector of dates for the data
    dates <- seq(from = as.Date("2021/1/1"),to = as.Date("2021/12/31") ,by = "day")
    
    # Create data containing simulated cases for each date and delay
    data <- data.frame(t = rep(1:365,times = 7),d = rep(1:7,each = 365),
                       Cases = as.vector(simcases),Date = rep(dates,times = 7))
    
    # Create a vector for the reporting date
    data <- data %>% mutate(Rep.date = Date+d)
    
    # Match data structure with vanDeKassteele et al. method
    epid1 <- data[,c("Date","Rep.date","Cases")]
    epi.data <- epid1 %>% 
      type.convert(as.is = TRUE) %>% 
      uncount(Cases)
    colnames(epi.data) <- c("onset.date", "report.date")
    epi.data$onset.date <- as.Date(epi.data$onset.date)
    epi.data$report.date <- as.Date(epi.data$report.date)
    
    # Generate prior reporting delay distribution
    f.priordelay <- genPriorDelayDist(mean.delay = 2, max.delay = 7, p = 0.99)
    
    ###### Data setup
    # Data setup
    rep.data <- dataSetup(
      data         = epi.data,
      start.date   = start.date, # Starting date of outbreak
      end.date     = NULL, # Ending date of outbreak (in real-time, leave NULL so end.date = nowcast.date)
      nowcast.date = nowcast.date, # Nowcast date
      days.back    = nowcast.date-start.date+1, # Number of days back from nowcast.date to include in estimation procedure
      f.priordelay = f.priordelay)          # Prior reporting delay PMF
    
    
    ###### Model setup
    # Model setup
    model.setup <- modelSetup(
      data = rep.data,
      ord = 2,
      kappa = list(u = 1e6, b = 1e6, s = 1e-6))
    
    ###### Nowcast
    # Nowcast
    nowcast.list <- nowcast(
      data = rep.data,
      model = model.setup,
      conf.level = 0.95)
    
    nowcast_result <- nowcast.list$nowcast %>% 
      filter(Date >= max(Date)-max.delay+1)
    obs <- rep.data %>% group_by(Date) %>%
      summarise(Obs=sum(Cases)) %>% 
      filter(Date>=max(Date)-max.delay+1) %>% 
      select(Obs)
    nowcast_result$Obs <- obs$Obs
    
    T.now <- as.numeric(nowcast.date - start.date) + 1
    dates.now <- seq(T.now-(max.delay-1),T.now) # dates that have nowcast
    mu.now <- f(a1  = a1, a2 = a2, x = dates.now)
    
    result <- nowcast_result %>%
      mutate(x = ifelse(Obs >= lwr & Obs <= upr,1,0)) %>%
      mutate(PI_width = upr - lwr,
             Bias_mu = med - mu.now,
             Relbias_mu = abs((med - mu.now) / mu.now))
    
    PI_cov[,j] <- result$x
    PI_width[,j] <- result$PI_width
    Bias_mu[,j] <- result$Bias_mu
    Relbias_mu[,j] <- result$Relbias_mu
    
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  print(j)
}


######## (a) PI COVERAGE #######
round(apply(PI_cov, MARGIN = 1, mean,na.rm = T)*100,2)

######## (b) PI width #######
round(apply(PI_width, MARGIN = 1, mean,na.rm = T),2)

######## (c) Bias(mu) #######
round(apply(Bias_mu, MARGIN = 1, mean,na.rm = T),2)

######## (d) Relative bias (mu) #######
round(apply(Relbias_mu, MARGIN = 1, mean,na.rm = T)*100,2)

results <- data.frame(rbind("PI Coverage" = round(apply(PI_cov, MARGIN = 1, mean, na.rm = TRUE) * 100, 2),
                            "PI Width" = round(apply(PI_width, MARGIN = 1, mean, na.rm = TRUE), 2),
                            "Bias (mu)" = round(apply(Bias_mu, MARGIN = 1, mean, na.rm = TRUE), 2),
                            "Relative bias (mu)" = round(apply(Relbias_mu, MARGIN = 1, mean, na.rm = TRUE) * 100, 2)))
colnames(results) <- c("T-6", "T-5", "T-4", "T-3", "T-2", "T-1", "T")

# Performance measures on the nowcast day
results[, 7, drop = FALSE]
