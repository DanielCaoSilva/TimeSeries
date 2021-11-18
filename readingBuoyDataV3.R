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
# Mauricio Missing data points ----
buoy_gaps <- buoy_ts %>%
  count_gaps(.full = TRUE)
buoy_ts <- buoy_ts %>% fill_gaps()

#Try fixing missing data using interpolate from ARIMA model---
buoy_fill <- buoy_ts %>% 
  model(ARIMA(V9 ~ pdq(4,1,0) + PDQ(2,0,0), approximation = FALSE)) %>%
  interpolate(buoy_ts)

# Box - Cox Transformation and Test Michael ----


#Unit Root Test - KSS Yodd ----


#Test Plot the data ----
buoy_ts %>% autoplot()
buoy_ts %>% gg_tsdisplay(V9,plot_type="partial",lag_max = 100)
buoy_ts %>% gg_tsdisplay(difference(V9),plot_type="partial",lag_max = 17)
buoy_ts %>% gg_tsdisplay(difference(V9,16),plot_type="partial",lag_max = 100)


#STL decomp (Needs work) Daniel ----
buoy_fill %>%
  model(
    STL(V9 ~ trend(window = 7)+
          season(window = "periodic"),
        robust = TRUE)) %>%
  components()%>%
  autoplot()
  

# Model Fitting ----
fit <- buoy_fill %>% model(auto = ARIMA(V9),
                          m1 = ARIMA(V9 ~ 0 + pdq(5,1,0) +PDQ(2,0,0,period = 24)),
                          m2 = ARIMA(V9 ~ 0 + pdq(4,1,0) +PDQ(2,0,0,period = 2)),
                          m3 = ARIMA(V9 ~ 0 + pdq(2,0,0) +PDQ(0,1,5,period = 16),
                                     order_constraint = p+q+P+Q<=10&(constant+d+D<=3)),
                          m4 = ARIMA(V9 ~ 0 + pdq() + PDQ(,period = 16))
                          )

# Residual Analysis----
fit %>% select(m3) %>% gg_tsresiduals(lag = 50)

augment(fit) %>%
  filter(.model == "m3") %>% 
  features(.innov,ljung_box,lag=50,dof=7)


# Forecast ----
fit %>% forecast(h="2 years") %>% autoplot(buoy_fill)
glance(fit)



########### Saved report(fit)
# Series: V9 
# Model: ARIMA(4,1,0)(2,0,0)[2] 
# 
# Coefficients:
#   ar1      ar2     ar3      ar4    sar1    sar2
# -0.2721  -0.1224  0.0637  -0.0759  0.1192  0.1390
# s.e.   0.0075   0.0327  0.0125   0.0173  0.0325  0.0195
# 
# sigma^2 estimated as 0.003467:  log likelihood=24831.08
# AIC=-49648.16   AICc=-49648.15   BIC=-49593.74







