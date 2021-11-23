# Daniel Silva
# Time Series Analysis on Significant Wave Height in the Santa Monica Bay
# Modeling with SARIMA and Dynamic Harmonic Regression

# Import libraries ----
library(astsa)
library(fpp3)
library(readr)
library(lubridate)

# Real time data (last 45 days)
#b_data <- read.table("https://www.ndbc.noaa.gov/data/realtime2/46268.txt")
#b_data <- read.table("https://www.ndbc.noaa.gov/data/realtime2/46268.txt")

# Load Data ----
#V1   2  3  4  5   6    7   8     9     10    11  12     13   14    15    16    17  18    19
#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS PTDY  TIDE
#yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC  nmi  hPa    ft
b_data <- read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=46221h2020.txt.gz&dir=data/historical/stdmet/")
b_data19 <- read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=46221h2019.txt.gz&dir=data/historical/stdmet/")



#Format the dates using a Date Class
b_date <- paste(b_data[,1],"-",b_data[,2],"-",b_data[,3]," ",b_data[,4],":",b_data[,5]) 
b_date <- ymd_hm(b_date,tz = "UTC")
b_date19 <- paste(b_data19[,1],"-",b_data19[,2],"-",b_data19[,3]," ",b_data19[,4],":",b_data19[,5]) 
b_date19 <- ymd_hm(b_date19,tz = "UTC")

#7693 to 7694 is when the change of intervals happens 17:26 to 19:53
#15416 to 15417 is when the change of intervals happens 20:23 to 22:26
b_data %>% select(V9) -> data
data$time <- b_date
b_data19 %>% select(V9) -> data19
data19$time <- b_date19



# Attempt to fix the uneven time intervals as select times ----
testies<-data$time[7694:15416]
fixedTimes <- testies + minutes(3)
data$time[7694:15416] <- fixedTimes

testies19<-data19$time[1:9906]
fixedTimes19 <- testies19 - minutes(4)
data19$time[1:9906] <- fixedTimes19

which((data$time[2:16184]-data$time[1:16183])>30)

newData <- rbind(data19,data)



# Turn into a tsibble ----
# using build_tsibble instead of as_tsibble to specify the interval  
buoy_ts <- build_tsibble(
  newData,
  key = dplyr::starts_with("2019-01-01 00:26:00"),
  index = time,
  ordered = TRUE,
  interval = new_interval(minute = 30),
  validate = TRUE
)




# Address additional missing data points (Maurico) ----
# TODO




# Count all the missing data and make them NA values ----
buoy_gaps <- buoy_ts %>%
  count_gaps(.full = TRUE)
buoy_ts <- buoy_ts %>% fill_gaps()



#Try fixing missing data using interpolate from ARIMA model----
buoy_fill <- buoy_ts %>% 
  model(ARIMA(V9 ~ pdq(4,1,0) + PDQ(2,0,0), approximation = FALSE)) %>%
  interpolate(buoy_ts)



# Box - Cox Transformation and Test (Michael) ----
# TODO



# Unit Root Test - KSS (Yodd) ----
# TODO



# Test Plot the data ----
buoy_ts %>% autoplot()
buoy_ts %>% gg_tsdisplay(V9,plot_type="partial",lag_max = 100)
buoy_ts %>% gg_tsdisplay(difference(V9),plot_type="partial",lag_max = 17)
buoy_ts %>% gg_tsdisplay(difference(V9,16),plot_type="partial",lag_max = 100)




#STL decomp ----
# Need to determine seaonal period(s)
lun_year <- 39000/30
lun_month <- (12*60+44+29*24*60)/30
lun_day <- (24+50/60)*2
sun_day <- 48
common_periods(buoy_fill)
buoy_fill %>%
  model(
    STL(log(V9) ~ season(period = lun_year)+       # 1 year
          season(period = lun_month)+  # months
          #season(period = 48*7)+   # weekly
          season(period = lun_day)+    # lunar daily
          season(period = sun_day),     # normal daily
        robust = TRUE)
  ) %>%
  components() %>%
  autoplot() + labs(x = "Observation")




# Model Fitting ----
# SARIMA model fit ----
fit <- buoy_fill %>% model(auto = ARIMA(V9),
                           m1 = ARIMA(V9 ~ 0 + pdq(5,1,0) +PDQ(2,0,0,period = 24)),
                           m2 = ARIMA(V9 ~ 0 + pdq(4,1,0) +PDQ(2,0,0,period = 2)),
                           m3 = ARIMA(V9 ~ 0 + pdq(2,0,0) +PDQ(0,1,5,period = 16),
                                      order_constraint = p+q+P+Q<=10&(constant+d+D<=3)),
                           m4 = ARIMA(V9 ~ 0 + pdq() + PDQ(,period = 16)),
)
# Dynamic Harmonic Regression fit ----
# Combines: Fourier terms for capturing seasonality with ARIMA errors capturing other dynamics 
# Need to know what the multiple seasonal periods are in order to fit
fitDHR <- model(buoy_fill,
                mk1 = ARIMA(log(V9) ~ PDQ(0,0,0)+pdq(d=0)+
                              fourier(period = lun_year, K = 3)+
                              fourier(period = lun_month, K = 1)+
                              fourier(period = lun_day, K = 1)+
                              fourier(period = sun_day, K = 2)
                ),
)

# Function for selection the optimal values of K for DHR+ARIMA model ----
kSelectors <- function(data_ts,
                       period1,
                       period2,
                       period3,
                       period4,
                       d_value,
                       k_start,
                       k_range){
  p1 <- period1
  p2 <- period2
  p3 <- period3
  p4 <- period4
  #k_list <- list(k1s = rep(NA), k2s = rep(NA), k3s = rep(NA), k4s = rep(NA))
  #aiccs <- list(rep(NA,k_range),rep(NA,k_range),rep(NA,k_range),rep(NA,k_range))
  tempSmallest <- 9999999
  smallestKvals <- rep(NA,4)
  modelIndex <- 0
  for (k1 in k_start:k_range) {
    print(k1)
    for (k2 in k_start:k_range) {
      print(k2)
      for (k3 in k_start:k_range) {
        for (k4 in k_start:k_range) {
          # Tests all the possible models using specific periods for DHR
          fits <- data_ts %>% 
            model(
              m1 = ARIMA(log(V9) ~ PDQ(0,0,0) + pdq(d=d_value) + #lun_year/month
                           fourier(period = p1, K = k1) + 
                           fourier(period = p2, K = k2)),
              m2 = ARIMA(log(V9) ~ PDQ(0,0,0) + pdq(d=d_value) + #lun_year/month/day
                           fourier(period = p1, K = k1) + 
                           fourier(period = p2, K = k2) + 
                           fourier(period = p3, K = k3)),
              m3 = ARIMA(log(V9) ~ PDQ(0,0,0) + pdq(d=d_value) + #lun_year/month + sun_day
                           fourier(period = p1, K = k1) + 
                           fourier(period = p2, K = k2) + 
                           fourier(period = p4, K = k4)),
              m4 = ARIMA(log(V9) ~ PDQ(0,0,0) + pdq(d=d_value) + #lun_year/month/day + sun_day
                           fourier(period = p1, K = k1) + 
                           fourier(period = p2, K = k2) + 
                           fourier(period = p3, K = k3) +
                           fourier(period = p4, K = k4)))
          
          # Takes the smallest AICc values of the associated model and saves it
          checkSmallest<-min(glance(fits)$AICc)
          whichModel <- which.min(glance(fits)$AICc)
          if(checkSmallest < tempSmallest){
            tempSmallest <- checkSmallest   # smallest AICc value
            smallestKvals <- c(k1,k2,k3,k4) # k values associated with smallest aicc
            modelIndex <- whichModel        # Which model gave the smallest AICc from the fit
          }
          #####
        }
      }
    }
  }
  # Return values in a list
  return_list <- list(smallestKvals,modelIndex,tempSmallest)
  return(return_list)
}

# Run the kSelectors function and assign the best fit to a new model to do analysis
bestKvals <- kSelectors(buoy_fill,
           period1 = lun_year,
           period2 = lun_month,
           period3 = lun_day,
           period4 = sun_day,
           d_value = 0,
           k_start = 1,
           k_range = 3)

#bestKvals[[1]]

fitDHR2 <- model(buoy_fill,
                mk1 = ARIMA(log(V9) ~ PDQ(0,0,0)+pdq(d=0)+
                              fourier(period = lun_year, K = 3)+
                              fourier(period = lun_month, K = 1)+
                              fourier(period = lun_day, K = 1)+
                              fourier(period = sun_day, K = 2)
                ),
)




# Residual Analysis ----
fit %>% select(m3) %>% gg_tsresiduals(lag = 50)
fitDHR %>% gg_tsresiduals()

augment(fit) %>%
  filter(.model == "m3") %>% 
  features(.innov,ljung_box,lag=50,dof=7)
augment(fitDHR) %>%
  features(.innov,ljung_box,lag=50)



# Training with Cross Validation ----
# TODO




# Forecast ----
fit %>% forecast(h = "2 years") %>% autoplot(buoy_fill)
glance(fit)
fitDHR %>% forecast(h = "3 months") %>% autoplot(buoy_fill)
glance(fitDHR)



# Error Analysis ----
# TODO








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




