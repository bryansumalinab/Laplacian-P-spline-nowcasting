library(blapsr)
library(Matrix)
library(numDeriv)
library(tidyverse)
library(RColorBrewer)

source("Nowcasting.R")

######### Data setup #########

data<-readxl::read_excel("mort2021.xlsx")

#day of the week (Monday as reference category)
data$Day=weekdays(as.Date(data$Rep.date))
data$Day <- relevel(factor(data$Day),"maandag")

date.now=as.Date('2021-10-31') #set nowcast date
max.delay=7 #maximum delay
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


########## Nowcasting #########

my_nowcast<-Nowcasting(data=data,date.now=date.now)


############# Plot nowcast  #############
plot_nowcast=function(data,date.now,data_CI){
  
  # Set colors
  blue <- brewer.pal(n = 9, name = "Set1")[2]
  oran <- brewer.pal(n = 9, name = "Set1")[5]
  grey <- brewer.pal(n = 9, name = "Set1")[9]
  
  D=max(data$d)
  data1 <- data %>%
    mutate(Reported = Reported %>%
             factor(
               levels = c("Reported", "Not yet reported", "Nowcast"),
               labels = c("Reported", "Not yet reported", "Nowcast"))) %>% 
    group_by(Date, Reported) %>%
    summarize(Cases = sum(Cases)) %>% 
    mutate(Date=as.Date(Date)) %>% as.data.frame() %>%
    filter(Date>date.now-(2*D))
  
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
    geom_crossbar(data=data_CI,
                  mapping=aes(x=Date,y=Cases,ymin=lower,ymax=upper),
                  fill = adjustcolor(oran, alpha = 0.3), colour = NA, width = 1,
                  inherit.aes = FALSE)
  
}

plot_nowcast(data=data,date.now=date.now,data_CI=my_nowcast$data_CI)


######## Delay density #########
my_nowcast$delay %>%
  ggplot(aes(Date , Delay, z = density))+
  geom_contour()+
  geom_contour_filled()+
  scale_x_date(date_breaks = "1 month",
               date_labels = "%b")
