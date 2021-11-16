library(astsa)
library(fpp3)
library(readr)
library(lubridate)

#wav_data <-read.csv("https://www.ndbc.noaa.gov/data/realtime2/46268.txt")
#b_data <- read.table("https://www.ndbc.noaa.gov/data/realtime2/46268.txt")
#b_data <- read.table("https://www.ndbc.noaa.gov/data/realtime2/46268.txt")
b_data <- read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=46221h2020.txt.gz&dir=data/historical/stdmet/")
#V1   2  3  4  5   6    7   8     9     10    11  12     13   14    15    16    17  18    19
#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS PTDY  TIDE
#yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC  nmi  hPa    ft


#Format the dates using a Date Class
b_date <- paste(b_data[,1],"-",b_data[,2],"-",b_data[,3]," ",b_data[,4],":",b_data[,5]) 
b_date <- ymd_hm(b_date,tz = "UTC")

#7693 to 7694 is when the change of intervals happens 17:26 to 19:53
#15416 to 15417 is when the change of intervals happens 20:23 to 22:26
b_date <- b_date
head(b_date)

b_data %>% select(V9) -> data
data$time <- b_date

#Attempt to fix the uneven time intervals as select times
testies<-data$time[7694:15416]
fixedTimes <- testies + minutes(3)
data$time[7694:15416] <- fixedTimes


which((data$time[2:16184]-data$time[1:16183])>30)

# Turn into a tsibble 
# using build_tsibble instead of as_tsibble to specifiy the interval  
buoy_ts <- build_tsibble(
  data,
  key = dplyr::starts_with("2020-01-01 00:26:00"),
  index = time,
  ordered = TRUE,
  interval = new_interval(minute = 30),
  validate = TRUE
)

# Count all the missing data and make them NA values
buoy_gaps <- buoy_ts %>%
  count_gaps(.full = TRUE)
buoy_ts <- buoy_ts %>% fill_gaps()

#Try fixing missing data using interpolate
buoy_fill <- buoy_ts %>% 
  model(ARIMA(V9)) %>%
  interpolate(buoy_ts)

#Test Plot the data 
buoy_ts %>% autoplot()
buoy_ts %>% gg_tsdisplay(V9,plot_type="partial")

#STL decomp 
buoy_fill %>%
  model(
    STL(V9 ~ trend(window = 7)+
          season(window = "periodic"),
        robust = TRUE)) %>%
  components()%>%
  autoplot()
  

# Model Fitting 
fit <- buoy_fill %>% model(auto = ARIMA(V9))
#                          m1 = ARIMA(log(V9) ~ 0 + pdq(1,1,0) +PDQ(2,0,2))+season("2") )
report(fit)
glance(fit)


# Forecast 
fit %>% forecast(h=1000) %>% autoplot(buoy_fill)
report(fit)










