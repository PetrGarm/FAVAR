library("readxl")
library("tidyverse")
library("forecast")
library("vars")
library("lmtest")
library("readxl")

setwd("C:/Users/petrg/Desktop/Диплом/Code/my experiments/")
source("my_tsCV.R")
#source("FAVAR_CV.R")

for_modeling <- read_excel("C:/Users/petrg/Desktop/Данные для диплома/real_big_data.xlsx", 
                            range = "R2C3:R242C296")
N_for_modeling <- dim(for_modeling)[2]

# Complete series only
complete_series_indices = c()
for(i in 1:N_for_modeling) {
  if (sum(is.na(for_modeling[,i])) == 0) {
    complete_series_indices = c(complete_series_indices, i)
  }
}

df <- for_modeling[,complete_series_indices]
colnames(df) <- colnames(for_modeling)[complete_series_indices]
ts <- ts(df, start=c(2000,1), frequency = 12)
N <- dim(ts)[2]
T <- dim(ts)[1]

transformations <- rep(0, N)

# Standardize & SA for some reason (correct PCA extraction probably)
for (i in 1:N){
  ts[,i] <- (ts[,i] - mean(ts[,i])) / sd(ts[,i])
  #[,i] <- ts[,i] - stl(ts[,i], 'periodic')$time.series[,1]
  transformations[i] <- ndiffs(ts[,i])
}

# induce stationarity
for (i in 1:N) {
  NAs <- rep(NA, transformations[i])
  if(transformations[i] != 0) {
    ts[,i] = c(NAs, diff(ts[,i], transformations[i]))
  }
}
max_I_order <- max(transformations)
ts <- na.omit(ts)

slow_index <- seq(N)[endsWith(colnames(ts), "_slow")]
y_names <- c("aRUCPI_slow", "aRUIP_slow")
y_index <- seq(N)[colnames(ts) %in% y_names]

'%not_in%' <- Negate('%in%')
x_index <- which(seq(N) %not_in% y_index)
slow_index_not_y <- slow_index[slow_index %not_in% y_index] # w.r.t ts object

Y <- ts[,y_index]
X <- ts[,x_index]
X_slow <- ts[,slow_index_not_y]
print(colnames(Y))
T <- dim(Y)[1]

# gdp
errors_FAVAR_cv <- my_tsCV(y = Y[,2], forecastfunction = fore_FAVAR, y_name="aRUIP_slow",
                           h = 12, X = X, Y = Y, X_slow = X_slow, K = 3, i_name = "RUCBIR=ECI", 
                           initial = 200)
res_FAVAR <- sqrt(colMeans(errors_FAVAR_cv^2, na.rm = TRUE))

### auto.ARIMA model CV
fore_arima <- function(y, h) {
  print(paste0(240 - length(y), " iterations left"))
  return(forecast(auto.arima(y) ,h = h))
}
errors_auto_arima_cv <- tsCV(Y[,3], forecastfunction = fore_arima, h = 12,
                             initial = 200)
res_ARIMA <- sqrt(colMeans(errors_auto_arima_cv^2, na.rm = TRUE))

# Auto ets model CV
fore_ets <- function(y, h) {
  print(paste0(240 - length(y), " iterations left"))
  return(forecast(ets(y) ,h = h))
}
errors_auto_ets_cv <- tsCV(Y[,3], forecastfunction = fore_ets, h = 12,
                           initial = 200)
res_ETS <- sqrt(colMeans(errors_auto_ets_cv^2, na.rm = TRUE))


# RW model CV
fore_rwf <- function(y, h) {
  model <- Arima(y,order = c(0,1,0), include.drift = TRUE)
  return(forecast(model, h = h))
}
errors_rw <- tsCV(Y[,3], forecastfunction = fore_rwf, h = 12,
                  initial = 200)
res_rw <- sqrt(colMeans(errors_rw^2, na.rm = TRUE))

results <- list(FAVAR = res_FAVAR, ARIMA = res_ARIMA, ETS = res_ETS, RWD = res_rw)
results <- t(data.frame(results))
View(results)
