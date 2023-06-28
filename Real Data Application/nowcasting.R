library(Matrix)
library(numDeriv)
library(dplyr)
library(ggplot2)


# Source functions
list.files(path = "LPS_functions", full.names = TRUE) %>% 
  purrr::walk(source)


######################################################################
######################### Mortality 2021 data ########################
data <- readxl::read_excel("mort2021.xlsx")
date.now <- as.Date('2021-03-31') # Nowcast date
max.delay <- 7

##### Nowcasting Poisson
nowcast.mortality.pois <- Nowcasting.Pois(data = data, # Data containing date of occurence and reporting
                                   date.now = date.now, # Nowcast date
                                   max.delay = max.delay, # Maximum delay
                                   day.effect = T, # True if include day of the week effect
                                   ref.day = "Monday" # reference category for day of the week effect
)

# Nowcast plot
plot.nowcast(nowcast.mortality.pois)
# Delay distribution plot
plot.delay(nowcast.mortality.pois)


##### Nowcasting NB
nowcast.mortality.NB <- Nowcasting.NB(data = data, # Data containing date of occurence and reporting
                                          date.now = date.now, # Nowcast date
                                          max.delay = max.delay, # Maximum delay
                                          day.effect = T, # True if include day of the week effect
                                          ref.day = "Monday" # reference category for day of the week effect
)

# Nowcast plot
plot.nowcast(nowcast.mortality.NB)
# Delay distribution plot
plot.delay(nowcast.mortality.NB)



######################################################################
######################### Incidence 2022 data ########################
data <- readxl::read_excel("incidence2022.xlsx")
date.now <- as.Date('2022-03-31') # Nowcast date
max.delay <- 5 # Maximum delay

##### Nowcasting Poisson
nowcast.mortality.pois <- Nowcasting.Pois(data = data, # Data containing date of occurence and reporting
                                          date.now = date.now, # Nowcast date
                                          max.delay = max.delay, # Maximum delay
                                          day.effect = T, # True if include day of the week effect
                                          ref.day = "Monday" # reference category for day of the week effect
)

# Nowcast plot
plot.nowcast(nowcast.mortality.pois)
# Delay distribution plot
plot.delay(nowcast.mortality.pois)


##### Nowcasting NB
nowcast.mortality.NB <- Nowcasting.NB(data = data, # Data containing date of occurence and reporting
                                      date.now = date.now, # Nowcast date
                                      max.delay = max.delay, # Maximum delay
                                      day.effect = T, # True if include day of the week effect
                                      ref.day = "Monday" # reference category for day of the week effect
)

# Nowcast plot
plot.nowcast(nowcast.mortality.NB)
# Delay distribution plot
plot.delay(nowcast.mortality.NB)

