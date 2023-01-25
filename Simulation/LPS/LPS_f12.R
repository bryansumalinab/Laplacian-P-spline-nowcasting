library(blapsr)
library(Matrix)
library(numDeriv)
library(tidyverse)
library(RColorBrewer)


source("Nowcasting_sim.R")

# Function f12
f <- function(a1,a2,x) exp(a1 + a2*sin((2*pi*x)/150))
mu.sim <- f(a1  = 3, a2 = 1, x = 1:365)

# Different nowcast dates (every end of the month) from March to November
dates.nowcast <- as.Date(c('2021-03-31','2021-04-30','2021-05-31','2021-06-30',
                           '2021-07-31','2021-08-31','2021-09-30','2021-10-31','2021-11-30'))

# Set the nowcast date
# Example for March 31, 2021
date.now <- dates.nowcast[1]

N <- 1000 # Number of simulations

p <- c(0.1,0.4,0.2,0.1,0.1,0.05,0.05) # Delay probabilities
max.delay <- 7 # Maximum delay

CI <- matrix(nrow = max.delay,ncol = N)
bias <- matrix(nrow = max.delay,ncol = N)
mape <- matrix(nrow = max.delay,ncol = N)
smape <- matrix(nrow = max.delay,ncol = N)

set.seed(12345)
for (j in 1:N) {
  tryCatch({
    
    simcases <- matrix(nrow = 365,ncol = 7)
    for (i in 1:365) {
      # Generate multinomial samples as explained in the simulation procedure
      simcases[i,] <- rmultinom(n = 1,size = rpois(n = 1,mu.sim[i]), prob = p)
    }
    
    # Create a vector of dates for the data
    dates <- seq(from = as.Date("2021/1/1"),to = as.Date("2021/12/31") ,by = "day")
    
    # create data containing simulated cases for each date and delay
    data <- data.frame(t = rep(1:365,times = 7),d = rep(1:7,each = 365),
                       Cases = as.vector(simcases),Date = rep(dates,times = 7))
    # Create a vector for the reporting date
    data <- data %>% mutate(Rep.date = Date+d)
    
    # Create numeric value for the nowcast day
    T.now <- as.numeric(date.now - as.Date(min(data$Date))) + 1
    
    # Distinguish types of cases (Reported, Not-yet-reported, Future cases)
    data <- data %>% mutate(Reported = ifelse(t + d <= T.now, "Reported",
                                              ifelse(t > T.now,"Future",
                                                     "Not yet reported")) %>%
                              factor(levels = c("Reported", "Not yet reported", "Future")))
    
    data <- data %>%
      filter(Reported != "Future") %>% # Exclude future cases
      mutate(Date = as.Date(Date), # Date the case occurred
             Rep.date = as.Date(Rep.date), # Reporting date
             Reported = Reported %>%
               factor(levels = c("Reported", "Not yet reported")
                      ,labels = c("Reported", "Not yet reported"))) %>%
      data.frame()
    
    # Do the nowcast using the Nowcasting_sim function
    my_nowcast <- Nowcasting_sim(data = data,date.now = date.now)
    
    CI[,j] <- my_nowcast$result$x
    bias[,j] <- my_nowcast$result$bias
    mape[,j] <- my_nowcast$result$mape
    smape[,j] <- my_nowcast$result$smape
    
  },error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # Continue with the loop when there is an error message
  # This happens when there are unusual simulated case
  print(j)
}



######## (a) CI COVERAGE #######
# Compute the prediction interval coverage for 1000 simulations
COVERAGE<- round(apply(CI, MARGIN = 1, mean,na.rm = T)*100,2)
COVERAGE

######## (b) MAPE #######
# Replace infinity with NA
# When there is zero simulated (actual) cases, it will lead to infinity values
mape[sapply(mape, is.infinite)] <- NA

# Compute MAPE for N simulations
MAPE <- round(apply(mape, MARGIN = 1, mean,na.rm = T)*100,2)
MAPE

######## (c) SMAPE #######
# Compute SMAPE for N simulations
SMAPE <- round(apply(smape, MARGIN = 1, mean,na.rm = T)*100,2)
SMAPE
