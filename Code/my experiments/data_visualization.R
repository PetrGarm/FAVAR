library("readxl")
library("tidyverse")
library("forecast")
library("vars")
library("lmtest")

library("feasts")
library("tsibbledata")
library("dplyr")
library("ggplot2")
library("lubridate")

for_modeling <- read_excel("C:/Users/petrg/Desktop/Данные для диплома/real_big_data_alina.xlsx", 
                           range = "R2C3:R242C296")
N_for_modeling <- dim(for_modeling)[2]

df_gdp <- for_modeling[(71):240,]

# Complete series only
complete_series_indices = c()
for(i in 1:N_for_modeling) {
  if (sum(is.na(df_gdp[,i])) == 0) {
    complete_series_indices = c(complete_series_indices, i)
  }
}

df_gdp <- df_gdp[,complete_series_indices]
colnames(df_gdp) <- colnames(for_modeling)[complete_series_indices]
ts_gdp <- ts(df_gdp, start=c(2005,11), frequency = 12)
N_gdp <- dim(ts_gdp)[2]
T_gdp <- dim(ts_gdp)[1]

transformations_gdp <- rep(0, N_gdp)

# Standardize & SA for some reason (correct PCA extraction probably)
for (i in 1:N_gdp){
  ts_gdp[,i] <- (ts_gdp[,i] - mean(ts_gdp[,i])) / sd(ts_gdp[,i])
  transformations_gdp[i] <- ndiffs(ts_gdp[,i])
}

# induce stationarity
for (i in 1:N_gdp) {
  NAs <- rep(NA, transformations_gdp[i])
  if(transformations_gdp[i] != 0) {
    ts_gdp[,i] = c(NAs, diff(ts_gdp[,i], transformations_gdp[i]))
  }
}
ts_gdp <- na.omit(ts_gdp)

slow_index <- seq(N_gdp)[endsWith(colnames(ts_gdp), "_slow")]
y_names <- c("RUGDP=ECI_slow", "RUCBIR=ECI")
y_index <- seq(N_gdp)[colnames(ts_gdp) %in% y_names]

'%not_in%' <- Negate('%in%')
x_index <- which(seq(N_gdp) %not_in% y_index)
slow_index_not_y <- slow_index[slow_index %not_in% y_index] # w.r.t ts_gdp object

Y_gdp <- ts_gdp[,y_index]
X <- ts_gdp[,x_index]
X_slow <- ts_gdp[,slow_index_not_y]
print(colnames(Y_gdp))
T_gdp <- dim(Y_gdp)[1]


times <- yearmon(seq(2005 + 10/12, 2020 - 1/12, by = 1 / 12))
df <- as.data.frame(ts_gdp)

move = ''
SLOW ='slow'
russia_ts = data.frame(Month=as.Date(character()),
                    name=character(),
                    Moving=character(),
                    value=character()) 
for(i in 1:dim(df)[2]) {
  name = colnames(df)[i]
  if(grepl(SLOW, name)) {
    move = 'Slow'
  }
  else {
    move = "Fast"
  }
  df_temp = data.frame(Month=times,
                       name=rep(name, length(times)),
                       Moving=rep(move, length(times)),
                       value=df[,i]) 
  russia_ts = rbind(russia_ts, df_temp)
}

russia <- as_tsibble(russia_ts, key = c(name, Moving), index = Month)
#russia[,3] = as.character(russia[,3])

features <- russia %>% features(value, feat_stl)

require(ggrepel)
russia %>%
  features(value, feat_stl) %>%
  ggplot(aes(x = trend_strength, y = seasonal_strength_year, 
             color = Moving)) +
  stat_density_2d(aes(fill = Moving, alpha = ..level..), bins = 5, geom = "polygon")  + 
  coord_equal() + 
  xlim(c(0,1)) + ylim(c(0,1)) +
  #facet_wrap(vars(Moving), nrow = 1) + 
  geom_point() + geom_point(shape = 1,colour = "black") + 
  labs(x = "Trend strength", y = "Seasonal strength") 

#
source("my_tsCV.R")
errors_FAVAR_cv_gdp_k1 <- my_tsCV(y = y, forecastfunction = fore_FAVAR, y_name=y_name,
                               h = 12, X = X, Y = Y_gdp, X_slow = X_slow, K = 1, i_name = i_name, 
                               initial = 144)
res_FAVAR_gdp_k1 <- sqrt(colMeans(errors_FAVAR_cv_gdp_k1^2, na.rm = TRUE))


