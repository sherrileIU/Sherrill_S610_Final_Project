# Elizabeth Sherrill
# STAT-S610 Final Project

library(readr)
library(lubridate)
library(ggplot2)

station_name = 'G131'
end_time = 2019
tR = 0.1
alpha = 20

fit_GPS_data(station_name, end_time, tR, alpha)

plot_GNSS_data(station_name, tR, alpha)