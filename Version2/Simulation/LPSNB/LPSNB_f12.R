library(Matrix)
library(numDeriv)
library(tidyverse)
library(RColorBrewer)

# Source functions
list.files(path = "LPS_functions", full.names = TRUE) %>% 
  purrr::walk(source)


# Function f12
f <- function(a1,a2,x) 50 + exp(a1 + a2*sin((2*pi*x)/150))
mu.sim <- f(a1  = 3, a2 = 2, x = 1:365)

# Different nowcast dates (every end of the month) from March to November
dates.nowcast <- as.Date(c('2021-03-31','2021-04-30','2021-05-31','2021-06-30',
                           '2021-07-31','2021-08-31','2021-09-30','2021-10-31','2021-11-30'))

# Set the nowcast date
# Example for March 31, 2021
date.now <- dates.nowcast[1]

N <- 500 # Number of simulations
p <- c(0, 0.1,0.4,0.2,0.1,0.1,0.05,0.05) # Delay probabilities
max.delay <- 7 # Maximum delay
overdisp <- 10 # overdispersion parameter

CI <- matrix(nrow = max.delay,ncol = N)
bias <- matrix(nrow = max.delay,ncol = N)
mape <- matrix(nrow = max.delay,ncol = N)
smape <- matrix(nrow = max.delay,ncol = N)
ciwdth <- matrix(nrow = max.delay,ncol = N)

set.seed(12345)
for (j in 1:N) {
  tryCatch({
    
    simcases <- matrix(nrow = 365,ncol = max.delay+1)
    for (i in 1:365) {
      # Generate multinomial samples as explained in the simulation procedure
      simcases[i,] <- rmultinom(n = 1,size = MASS::rnegbin(n = 1, mu = mu.sim[i], theta = overdisp), prob = p)
    }
    
    # Create a vector of dates for the data
    dates <- seq(from = as.Date("2021/1/1"),to = as.Date("2021/12/31") ,by = "day")
    
    # Create data containing simulated cases for each date and delay
    data <- data.frame(t = rep(1:365,times = max.delay+1),d = rep(0:max.delay,each = 365),
                       Cases = as.vector(simcases),Date = rep(dates,times = max.delay+1))
    
    # Create a vector for the reporting date
    data <- data %>% mutate(Rep.date = Date+d)
    
    
    # Nowcasting
    nowcast <- Nowcasting.NB(data = data, # Data containing date of occurence and reporting
                             date.now = date.now, # Nowcast date
                             max.delay = max.delay, # Maximum delay
                             day.effect = F, # True if include day of the week effect
    )
    
    Obs <- nowcast$data %>% filter(Date %in% nowcast$nowcast$Date) %>% 
      group_by(Date) %>%
      summarise(ObsCases = sum(Cases))
    
    result <- nowcast$nowcast
    result$ObsCases <- Obs$ObsCases
    
    value <- result %>%
      mutate(x = ifelse(ObsCases >= `CI 95% lower` & ObsCases <= `CI 95% upper`,1,0)) %>%
      mutate(bias = Cases - ObsCases,
             mape = abs((Cases-ObsCases)/ObsCases),
             smape = abs(Cases-ObsCases)/(((abs(Cases)+abs(ObsCases)))/2),
             ciwdth = `CI 95% upper` - `CI 95% lower`)
    CI[,j] <- value$x
    bias[,j] <- value$bias
    mape[,j] <- value$mape
    smape[,j] <- value$smape
    ciwdth[,j] <- value$ciwdth
    
  },error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # Continue with the loop when there is an error message
  # This happens when there are unusual simulated case
  print(j)
}



######## (a) CI COVERAGE #######
# Prediction interval coverage for N simulations
COVERAGE<- round(apply(CI, MARGIN = 1, mean,na.rm = T)*100,2)
COVERAGE

######## (b) MAPE #######
# Replace infinity with NA
# When there is zero simulated (actual) cases, it will lead to infinity values
mape[sapply(mape, is.infinite)] <- NA

# Average MAPE for N simulations
MAPE <- round(apply(mape, MARGIN = 1, mean,na.rm = T)*100,2)
MAPE

######## (c) SMAPE #######
# Average SMAPE for N simulations
SMAPE <- round(apply(smape, MARGIN = 1, mean,na.rm = T)*100,2)
SMAPE

######## (D) CI width #######
# Average CI width for N simulations
CIwidth <- round(apply(ciwdth, MARGIN = 1, mean,na.rm = T),2)
CIwidth