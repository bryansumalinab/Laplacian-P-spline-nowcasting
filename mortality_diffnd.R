library(blapsr)
library(Matrix)
library(numDeriv)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)

Sys.setlocale("LC_TIME", "English") #set language to English

source("Nowcasting.R")

data <- readxl::read_excel("mort2021.xlsx")
data$Day <- weekdays(as.Date(data$Rep.date))
data$Day <- relevel(factor(data$Day),"Monday")

dates.now <- as.Date(c('2021-03-31','2021-04-30','2021-05-31','2021-06-30',
                    '2021-07-31','2021-08-31','2021-09-30','2021-10-31'))

max.delay <- 7 #maximum delay

f_data <- function(data,max.delay,date.now){
  date.start <- as.Date(min(data$Date))
  T.now <- as.numeric(date.now - date.start) + 1
  data <- data %>% mutate(Reported = ifelse(t + d <= T.now, "Reported",
                                         ifelse(t > T.now,"Future",
                                                "Not yet reported")) %>%
                         factor(levels = c("Reported", "Not yet reported", "Future")))
  
  data <- data %>% 
    filter(Reported != "Future") %>% 
    filter(d <= max.delay) %>%
    filter(Date <= date.now) %>% 
    mutate(Date = as.Date(Date),
           Rep.date = as.Date(Rep.date),
           Reported = Reported %>%
             factor(levels = c("Reported", "Not yet reported")
                    ,labels = c("Reported", "Not yet reported"))) %>% 
    data.frame()
  data
}


my_data <- list()
my_data <- lapply(dates.now,
                f_data,
                max.delay = max.delay,
                data = data)


my_nowcast <- mapply(Nowcasting,
                   my_data,
                   dates.now,
                   SIMPLIFY = F)

############# Plot nowcast  #############
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

myplots <- list()
for(i in 1:length(dates.now)){
  myplots[[i]] <- plot_nowcast(data = my_data[[i]],
                             data_CI = my_nowcast[[i]]$data_CI,
                             date.now = dates.now[i])
}

do.call("grid.arrange",c(myplots, ncol = 2))
# do.call("grid.arrange",c(myplots[1:4], ncol = 2))
# do.call("grid.arrange",c(myplots[5:8], ncol = 2))

######## Delay density #########
delay_plots <- list()

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

for(i in 1:length(dates.now)){
  delay_plots[[i]] <- plot.delay(delay = my_nowcast[[i]]$delay,
                               date.now = dates.now[i])
}

do.call("grid.arrange",c(delay_plots, ncol = 2))
# do.call("grid.arrange",c(delay_plots[1:4], ncol = 2))
# do.call("grid.arrange",c(delay_plots[5:8], ncol = 2))