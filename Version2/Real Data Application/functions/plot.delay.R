######## Delay density #########
plot.delay <- function(nowcast){
  nowcast$delay %>%
    ggplot(aes(Date , Delay, z = density))+
    geom_contour()+
    geom_contour_filled()+
    scale_x_date(date_breaks = "1 month",
                 date_labels = "%b")
}
  