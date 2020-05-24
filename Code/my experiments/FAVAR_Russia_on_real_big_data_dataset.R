library("readxl")
library("tidyverse")
library("forecast")
library("vars")
library("lmtest")

setwd("C:/Users/petrg/Desktop/Диплом/Code/my experiments/")
source("my_tsCV.R")
#source("FAVAR_CV.R")

for_modeling <- read_excel("C:/Users/petrg/Desktop/Данные для диплома/real_big_data_alina.xlsx", 
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
#lns <- rep(0, N)
# Standardize & log
for (i in 1:N){
  #if (means[i] > 10^4 & sum(ts[,i] <=0) == 0) {
    #lns[i] <- 1
    #ts[,i] <- log(ts[,i])
  #}
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
ts <- na.omit(ts)

slow_index <- seq(N)[endsWith(colnames(ts), "_slow")]
y_names <- c("aRUCPI_slow", "RUCBIR=ECI")
#y_names <- c("aRUCPI_slow", "RUCBIR=ECI", "RUUNR=ECI_slow")
y_index <- seq(N)[colnames(ts) %in% y_names]

'%not_in%' <- Negate('%in%')
x_index <- which(seq(N) %not_in% y_index)
slow_index_not_y <- slow_index[slow_index %not_in% y_index] # w.r.t ts object

Y <- ts[,y_index]
X <- ts[,x_index]
X_slow <- ts[,slow_index_not_y]
print(colnames(Y))
T <- dim(Y)[1]

# cpi
errors_FAVAR_cv <- my_tsCV(y = Y[,2], forecastfunction = fore_FAVAR, y_name="aRUCPI_slow",
                           h = 12, X = X, Y = Y, X_slow = X_slow, K = 1, i_name = "RUCBIR=ECI", 
                           initial = 200)
res_FAVAR <- sqrt(colMeans(errors_FAVAR_cv^2, na.rm = TRUE))


errors_VAR_cv <- my_tsCV(y = Y[,2], forecastfunction = fore_FAVAR, y_name="aRUCPI_slow",
                         h = 12, X = X, Y = Y, X_slow = X_slow, K = 1, i_name = "RUCBIR=ECI", use_VAR = T, 
                         initial = 200)
res_VAR <- sqrt(colMeans(errors_VAR_cv^2, na.rm = TRUE))

### auto.ARIMA model CV
fore_arima <- function(y, h) {
  print(paste0(240 - length(y), " iterations left"))
  return(forecast(auto.arima(y) ,h = h))
}
errors_auto_arima_cv <- tsCV(Y[,2], forecastfunction = fore_arima, h = 12,
                             initial = 200)
res_ARIMA <- sqrt(colMeans(errors_auto_arima_cv^2, na.rm = TRUE))

# Auto ets model CV
fore_ets <- function(y, h) {
  print(paste0(240 - length(y), " iterations left"))
  return(forecast(ets(y) ,h = h))
}
errors_auto_ets_cv <- tsCV(Y[,2], forecastfunction = fore_ets, h = 12,
                           initial = 200)
res_ETS <- sqrt(colMeans(errors_auto_ets_cv^2, na.rm = TRUE))


# RW model CV
fore_rwf <- function(y, h) {
  model <- Arima(y,order = c(0,1,0), include.drift = TRUE)
  return(forecast(model, h = h))
}
errors_rw <- tsCV(Y[,2], forecastfunction = fore_rwf, h = 12,
                  initial = 200)
res_rw <- sqrt(colMeans(errors_rw^2, na.rm = TRUE))


results <- list(FAVAR = res_FAVAR, ARIMA = res_ARIMA, ETS = res_ETS, VAR = res_VAR, RWD = res_rw)
results <- t(data.frame(results))
View(results)


### unemp
#y_names <- c("RUUNR=ECI_slow", "RUCBIR=ECI")
y_names <- c("aRUCPI_slow", "RUUNR=ECI_slow")
y_index <- seq(N)[colnames(ts) %in% y_names]

slow_index <- seq(N)[endsWith(colnames(ts), "_slow")]
'%not_in%' <- Negate('%in%')
x_index <- which(seq(N) %not_in% y_index)
slow_index_not_y <- slow_index[slow_index %not_in% y_index] # w.r.t ts object

Y <- ts[,y_index]
X <- ts[,x_index]
X_slow <- ts[,slow_index_not_y]
print(colnames(Y))
T <- dim(Y)[1]
y_name <- "RUUNR=ECI_slow"
i_name <- "aRUCPI_slow"
y <- Y[,2]

source("my_tsCV.R")
errors_FAVAR_cv_unemp <- my_tsCV(y = y, forecastfunction = fore_FAVAR, y_name=y_name,
                           h = 12, X = X, Y = Y, X_slow = X_slow, K = 1, i_name = i_name, 
                           initial = 200)
res_FAVAR_unemp <- sqrt(colMeans(errors_FAVAR_cv_unemp^2, na.rm = TRUE))

errors_VAR_cv_unemp <- my_tsCV(y = y, forecastfunction = fore_FAVAR, y_name=y_name,
                         h = 12, X = X, Y = Y, X_slow = X_slow, K = 1, i_name = i_name, use_VAR = T, 
                         initial = 200)
res_VAR_unemp <- sqrt(colMeans(errors_VAR_cv_unemp^2, na.rm = TRUE))

errors_auto_arima_cv_unemp <- tsCV(y, forecastfunction = fore_arima, h = 12,
                             initial = 200)
res_ARIMA_unemp <- sqrt(colMeans(errors_auto_arima_cv_unemp^2, na.rm = TRUE))

errors_auto_ets_cv_unemp <- tsCV(Y[,2], forecastfunction = fore_ets, h = 12,
                           initial = 200)
res_ETS_unemp <- sqrt(colMeans(errors_auto_ets_cv_unemp^2, na.rm = TRUE))

errors_rw_unemp <- tsCV(y, forecastfunction = fore_rwf, h = 12,
                  initial = 200)
res_rw_unemp <- sqrt(colMeans(errors_rw_unemp^2, na.rm = TRUE))

results_unemp <- list(FAVAR = res_FAVAR_unemp, ARIMA = res_ARIMA_unemp, ETS = res_ETS_unemp, 
                      VAR = res_VAR_unemp, RWD = res_rw_unemp)
results_unemp <- t(data.frame(results_unemp))
View(results_unemp)

#gdp
obs_gdp <- dim(na.omit(for_modeling['RUGDP=ECI_slow']))[1]
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
y_name <- "RUGDP=ECI_slow"
i_name <- "RUCBIR=ECI"
y <- Y_gdp[,2]


### forecast gdp
source("my_tsCV.R")
errors_FAVAR_cv_gdp <- my_tsCV(y = y, forecastfunction = fore_FAVAR, y_name=y_name,
                                 h = 12, X = X, Y = Y_gdp, X_slow = X_slow, K = 2, i_name = i_name, 
                                 initial = 144)
res_FAVAR_gdp <- sqrt(colMeans(errors_FAVAR_cv_gdp^2, na.rm = TRUE))

errors_VAR_cv_gdp <- my_tsCV(y = y, forecastfunction = fore_FAVAR, y_name=y_name,
                             h = 12, X = X, Y = Y_gdp, X_slow = X_slow, K = 1, i_name = i_name, 
                             use_VAR = T, initial = 144)
res_VAR_gdp <- sqrt(colMeans(errors_VAR_cv_gdp^2, na.rm = TRUE))

errors_auto_arima_cv_gdp <- tsCV(y, forecastfunction = fore_arima, h = 12,
                                   initial = 144)
res_ARIMA_gdp <- sqrt(colMeans(errors_auto_arima_cv_gdp^2, na.rm = TRUE))

errors_auto_ets_cv_gdp <- tsCV(y, forecastfunction = fore_ets, h = 12,
                                 initial = 144)
res_ETS_gdp <- sqrt(colMeans(errors_auto_ets_cv_gdp^2, na.rm = TRUE))

errors_rw_gdp <- tsCV(y, forecastfunction = fore_rwf, h = 12,
                        initial = 144)
res_rw_gdp <- sqrt(colMeans(errors_rw_gdp^2, na.rm = TRUE))

results_gdp <- list(FAVAR = res_FAVAR_gdp, ARIMA = res_ARIMA_gdp, ETS = res_ETS_gdp, 
                    VAR = res_VAR_gdp, RWD = res_rw_gdp)
results_gdp <- t(data.frame(results_gdp))
View(results_gdp)



### LATEX table generation
variables <- colnames(ts_gdp)
var_diffs <- transformations_gdp
data_info <- read_excel("C:/Users/petrg/Desktop/Данные для диплома/real_big_data_alina_info.xlsx", 
                        range = "A4:L298")

for (i in 1:length(variables)) {
  transform_code <- NA
  info <- sub("%","\\\\%",sub("Russia, "," ",data_info[data_info$RIC == sub("_slow","",variables[i]),]$Name))
  if(var_diffs[i] == 0)
    transform_code <- '$\\l$'
   if(var_diffs[i] == 1)
     transform_code <- paste0('$\\Delta','$')
   if(var_diffs[i] == 2) 
   transform_code <- paste0('$\\Delta^',var_diffs[i],'$')
  
  cat(paste0(i,". ",sub('_slow', '*',variables[i]), ' & ',info ,' & ', transform_code, '\\\\', '\n'))
}


library(xtable)
cpi_latex <- xtable(t(results*100), caption = 'RMSE $\\times$ 100 for CPI forecasts')
print(cpi_latex,hline.after=c(-1, 0), tabular.environment = "longtable")

unemp_latex <- xtable(t(results_unemp*100), caption = 'RMSE $\\times$ 100 for unemployment rate forecasts')
print(unemp_latex,hline.after=c(-1, 0), tabular.environment = "longtable")

gdp_latex <- xtable(t(results_gdp*100), caption = 'RMSE $\\times$ 100 for GDP forecasts')
print(gdp_latex,hline.after=c(-1, 0), tabular.environment = "longtable")
