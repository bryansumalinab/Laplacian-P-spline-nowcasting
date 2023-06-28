############# Plot nowcast  #############
plot.nowcast <- function(nowcast){
  data <- nowcast$data
  date.now <- nowcast$date.now
  data_CI <- nowcast$nowcast
  
  D <- max(data$d)
  data1 <- data %>%
    mutate(Reported = Reported %>%
             factor(
               levels = c("Reported", "Not yet reported", "Nowcast"),
               labels = c("Reported", "Not yet reported", "Nowcast"))) %>% 
    filter(Date>date.now-(2*D)) %>% 
    group_by(Date, Reported) %>%
    summarize(Cases = sum(Cases), .groups = "drop")
  
  
  ggplot(
    data = data1,
    mapping = aes(x = Date, y = Cases, fill = Reported)) +
    geom_col(
      width = 1,
      position = position_stack(reverse = TRUE)) +
    geom_step(
      data = data_CI,
      mapping = aes(x = Date + 0.5 , y = Cases),
      colour = "red",
      size = 0.8,
      direction = "vh",
      inherit.aes = FALSE) +
    scale_fill_manual(
      values = c("blue", "grey50", "red"),
      name = "",
      drop = FALSE)+
    geom_crossbar(data = data_CI,
                  mapping = aes(x = Date,y = Cases,ymin = `CI 95% lower`,ymax = `CI 95% upper`),
                  fill = adjustcolor("orange", alpha = 0.3), colour = NA, width = 1,
                  inherit.aes = FALSE)
}

