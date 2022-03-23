library(fable)
library(lubridate)
library(dplyr)
library(feasts)
library(ggplot2)
xm <-read.csv("Data/Donnees1.csv")
# Garder les dernières valeurs
xme <- ts(xm[1:(nrow(xm)-4),], start = 2000, frequency = 12) %>%
  as_tsibble()
xm_complete = ts(xm[,1], start = 2000, frequency = 12) %>%
  as_tsibble()

xme %>%
  gg_subseries()

model = xme  %>%
  model(arima1 = ARIMA((value) ~ 0 + pdq(0, 0, 1) + PDQ(0,1,0) ),
        arima2 = ARIMA((value) ~ 0 + pdq(1, 0, 1) + PDQ(0,1,0) ))
model %>% 
  coef()
model %>%
  accuracy()
p <- model %>%
  forecast(h="4 months") %>% 
  autoplot(xme  %>% 
              filter(year(index) >2010)) 
p
p +
  geom_line(data= tail(xm_complete, 5), aes(x = index, y = value),
            linetype = "dotted")

