library(blapsr)
library(Matrix)
library(numDeriv)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)

Sys.setlocale("LC_TIME", "English") # Set language to English

source("Nowcasting.R") # Load the Nowcasting.R script

data <- readxl::read_excel("mort2021.xlsx")

# Different nowcast dates
dates.now <- as.Date(c('2021-03-31','2021-04-30','2021-05-31','2021-06-30',
                       '2021-07-31','2021-08-31','2021-09-30','2021-10-31'))

# Day of the week (Monday as reference category)
data$Day <- weekdays(as.Date(data$Rep.date))
data$Day <- relevel(factor(data$Day),"Monday")

max.delay <- 7 # Maximum delay

# Function to create data that only include cases 
# up to the nowcast date and a specified maximum delay
f_data <- function(data, max.delay, date.now){
  T.now <- as.numeric(date.now - as.Date(min(data$Date))) + 1 # create numeric value for the nowcast day
  
  # Distinguish types of cases (Reported, Not-yet-reported, Future cases)
  data <- data %>% mutate(Reported = ifelse(t + d <= T.now, "Reported",
                                            ifelse(t > T.now,"Future",
                                                   "Not yet reported")) %>%
                            factor(levels = c("Reported", "Not yet reported", "Future")))
  
  data <- data %>% 
    filter(Reported != "Future") %>% 
    filter(d <= max.delay) %>%
    mutate(Date = as.Date(Date),
           Rep.date = as.Date(Rep.date),
           Reported = Reported %>%
             factor(levels = c("Reported", "Not yet reported")
                    ,labels = c("Reported", "Not yet reported"))) %>% 
    data.frame()
  data
}


my_data <- list()

# Apply the f_data() function for the different nowcast dates
my_data <- lapply(dates.now,
                  f_data,
                  max.delay = max.delay,
                  data = data)

# Do nowcasting for the different nowcast dates
my_nowcast <- mapply(Nowcasting,
                     my_data,
                     dates.now,
                     SIMPLIFY = F)

############# Plot nowcast  #############

# Function to provide plots for each nowcast dates
plot_nowcast <- function(data,date.now,data_CI){
  
  # Set colors
  blue <- brewer.pal(n = 9, name = "Set1")[2]
  oran <- brewer.pal(n = 9, name = "Set1")[5]
  grey <- brewer.pal(n = 9, name = "Set1")[9]
  
  D <- max(data$d)
  data1 <- data %>%
    mutate(Reported = Reported %>%
             factor(
               levels = c("Reported", "Not yet reported", "Nowcast"),
               labels = c("Reported", "Not yet reported", "Nowcast"))) %>% 
    group_by(Date, Reported) %>%
    summarize(Cases = sum(Cases), .groups = "drop") %>% 
    mutate(Date = as.Date(Date)) %>% as.data.frame() %>%
    filter(Date > date.now-(2*D))
  
  ggplot(
    data = data1,
    mapping = aes(x = Date, y = Cases, fill = Reported)) +
    geom_col(
      width = 1,
      position = position_stack(reverse = TRUE)) +
    geom_step(
      data = data_CI,
      mapping = aes(x = Date + 0.5 , y = Cases),
      colour = oran,
      size = 0.5,
      direction = "vh",
      inherit.aes = FALSE) +
    scale_fill_manual(
      values = c(blue, grey, oran),
      name = "",
      drop = FALSE)+
    geom_crossbar(data = data_CI,
                  mapping = aes(x = Date,y = Cases,ymin = lower,ymax = upper),
                  fill = adjustcolor(oran, alpha = 0.3), colour = NA, width = 1,
                  inherit.aes = FALSE)+
    theme(legend.position = "none") +
    labs(x = NULL, y = NULL)+
    ggtitle(paste("Nowcast date: ",toupper(format(date.now,"%b %d")) ))
  
}

# Create a list object for the nowcast plot of different
# nowcast date using the plot_nowcast() function
myplots <- list()
for(i in 1:length(dates.now)){
  myplots[[i]] <- plot_nowcast(data = my_data[[i]],
                               data_CI = my_nowcast[[i]]$data_CI,
                               date.now = dates.now[i])
}

# Plot the result
do.call("grid.arrange",c(myplots, ncol = 2))


######## Delay density #########

# Function to produce delay density plots for each nowcast dates
plot.delay <- function(delay,date.now){
  delay %>%
    ggplot(aes(Date , Delay, z = density))+
    geom_contour()+
    geom_contour_filled(breaks = seq(0,.8,by = 0.1))+
    scale_x_date(date_breaks = "1 month",
                 date_labels = "%b",
                 expand = c(0,0))+
    labs(x = NULL)+
    ggtitle(paste("Nowcast date: ",toupper(format(date.now,"%b %d")) ))
}

# Create a list object for the density plot of different
# nowcast date using the plot.delay() function
delay_plots <- list()
for(i in 1:length(dates.now)){
  delay_plots[[i]] <- plot.delay(delay = my_nowcast[[i]]$delay,
                                 date.now = dates.now[i])
}

# Plot the result
do.call("grid.arrange",c(delay_plots, ncol = 2))
