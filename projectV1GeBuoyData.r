library(astsa)
library(fpp3)
library(readr)
library(lubridate)

#wav_data <-read.csv("https://www.ndbc.noaa.gov/data/realtime2/46268.txt")
#b_data <- read.table("https://www.ndbc.noaa.gov/data/realtime2/46268.txt")
b_data <- read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=46221h2020.txt.gz&dir=data/historical/stdmet/")
#V1   2  3  4  5   6    7   8     9     10    11  12     13   14    15    16    17  18    19
#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS PTDY  TIDE
#yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC  nmi  hPa    ft

b_date = paste(b_data[,1],"-",b_data[,2],"-",b_data[,3]," ",b_data[,4],":",b_data[,5],":",00,sep="") 
#b_date = paste(b_data[,1],"-",b_data[,2],"-",b_data[,3]," ",b_data[,4],":",00,sep="") 

b_date = ymd_hms(b_date, tz = "UTC")
head(b_date)

b_data %>% select(V9) -> data
data$Time <- b_date

buoy_ts <-as_tsibble(data,regular = FALSE)

buoy_gaps <- buoy_ts %>%
  count_gaps(.full = TRUE)
# buoy_ts <- buoy_ts %>% fill_gaps(.full = TRUE)

buoy_ts %>% autoplot()

buoy_ts %>% mutate(log = log(V9)) %>%
  mutate(diff = difference(V9)) %>%
  mutate(diff2 = difference(diff)) %>%
  mutate(logDif = difference(log)) %>%
  mutate(logDif2 = difference(logDif)) -> buoy_ts
buoy_ts %>% gg_tsdisplay(V9,"partial")
buoy_ts %>% gg_tsdisplay(log,"partial")
buoy_ts %>% gg_tsdisplay(diff,"partial")
buoy_ts %>% gg_tsdisplay(diff2,"partial")
buoy_ts %>% gg_tsdisplay(logDif,"partial")
buoy_ts %>% gg_tsdisplay(logDif2,"partial")

# fit2 <-buoy_ts %>%
#   model(m1 = ARIMA(log ~ pdq(1,0,1)),
#         m2 = ARIMA(log ~ pdq(2,0,0)),
#         m3 = ARIMA(log ~ pdq(2,0,1)),
#         m4 = ARIMA(log ~ pdq(3,0,2)),)
# glance(fit2) %>% arrange(AICc)
fit <- buoy_ts %>% 
  model(arima = ARIMA(log ~ 0 + pdq(0,1,1) +PDQ(0,1,1))) %>% 
  report()

#wav_data_time %>% transform(test2,time=paste0(V1,V2,V3))  








