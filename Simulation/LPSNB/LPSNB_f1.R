rm(list = ls())

library(Matrix)

source('cubicbs.R')
source('nowcasting.NB.R')

# Function f1
a1 <- 3
a2 <- 1
f <- function(a1,a2,x) exp(a1 + a2*sin((2*pi*x)/150))
mu.sim <- f(a1  = a1, a2 = a2, x = 1:365)

# Different nowcast dates (every end of the month) from March to November
dates.nowcast <- as.Date(c('2021-03-31','2021-04-30','2021-05-31','2021-06-30',
                           '2021-07-31','2021-08-31','2021-09-30','2021-10-31','2021-11-30'))

# Set the nowcast date
# Example for March 31, 2021
date.now <- dates.nowcast[1]

p <- c(0,0.1,0.4,0.2,0.1,0.1,0.05,0.05) # delay probabilities
N <- 500 # number of simulations
max.delay <- 7 # maximum delay
overdisp <- 10

PI_cov <- matrix(nrow = max.delay,ncol = N)
PI_width <- matrix(nrow = max.delay,ncol = N)
Bias_mu <- matrix(nrow = max.delay,ncol = N)
Relbias_mu <- matrix(nrow = max.delay,ncol = N)

set.seed(12345)
for (j in 1:N) {
  tryCatch({
    
    simcases <- matrix(nrow = 365,ncol = max.delay+1)
    for (i in 1:365) {
      simcases[i,] <- rmultinom(n = 1,size = MASS::rnegbin(n = 1, mu = mu.sim[i], theta = overdisp), prob = p)
    }
    
    dates <- seq(from = as.Date("2021/1/1"),to = as.Date("2021/12/31") ,by = "day")
    
    data <- data.frame(t = rep(1:365,times = max.delay+1),d = rep(0:max.delay,each = 365),
                       Cases = as.vector(simcases),Date = rep(dates,times = max.delay+1))
    data$Rep.date <- data$Date + data$d
    
    
    # # Create a vector for the reporting date
    T.now <- as.numeric(date.now - as.Date(min(data$Date))) + 1
    data <- data[data$t <= T.now, ]
    Reported <- ifelse(data$t + data$d <= T.now, "Reported", "Not yet reported")
    Reported <- factor(Reported, levels = c("Reported", "Not yet reported"))
    data$Reported <- Reported
    
    # Nowcasting
    nowcast <- nowcasting.NB(data = data)
    
    # Filter data based on non-missing dates
    nowcast_dates <- nowcast$data[!is.na(match(nowcast$data$Date, na.omit(nowcast$cases.now)$Date)), ]
    # Group by Date and calculate the sum of Cases
    Obs <- aggregate(Cases ~ Date, data = nowcast_dates, FUN = sum)
    names(Obs)[2] <- "ObsCases"

    result <- na.omit(nowcast$cases.now)
    result$ObsCases <- Obs$ObsCases
    
    dates.now <- unique(nowcast_dates$t)
    mu.now <- f(a1  = a1, a2 = a2, x = dates.now)
    
    PI_cov[,j] <- ifelse(result$ObsCases >= result$CI95L & result$ObsCases <= result$CI95R, 1, 0)
    PI_width[,j] <- result$CI95R - result$CI95L
    Bias_mu[,j] <- result$y - mu.now
    Relbias_mu[,j] <- abs((result$y - mu.now) / mu.now)
  },error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
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