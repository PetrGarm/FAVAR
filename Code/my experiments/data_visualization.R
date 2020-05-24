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

source("FAVAR_CV.R")
source("my_tsCV.R")

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
y_name <- "RUGDP=ECI_slow"
i_name <- "RUCBIR=ECI"
y <- Y_gdp[,2]

### Visualization 199 ts ###
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

####

# gdp
errors_FAVAR_cv_gdp_k1 <- my_tsCV(y = y, forecastfunction = fore_FAVAR, y_name=y_name,
                               h = 12, X = X, Y = Y_gdp, X_slow = X_slow, K = 1, i_name = i_name, 
                               initial = 144)
res_FAVAR_gdp_k1 <- sqrt(colMeans(errors_FAVAR_cv_gdp_k1^2, na.rm = TRUE))

results_gdp = c()
for (k in 0:5) {
  errors_FAVAR_cv_gdp_k <- my_tsCV(y = y, forecastfunction = fore_FAVAR, y_name=y_name,
                                    h = 12, X = X, Y = Y_gdp, X_slow = X_slow, K = k + (k==0), 
                                    i_name = i_name, 
                                    initial = 144, use_VAR = (k==0))
  res_FAVAR_gdp_k <- sqrt(colMeans(errors_FAVAR_cv_gdp_k^2, na.rm = TRUE))
  results_gdp = c(results_gdp, sum(res_FAVAR_gdp_k))
}

df_gdp = data.frame(k=seq(0,5), sum_rmse=results_gdp)

ggplot(data=df_gdp, aes(x=k, y=sum_rmse)) + geom_point(shape=15, size=3) + geom_line() + labs(x='K', y='Sum of RMSE') +
  ggtitle("GDP forecast error, Russia")

# cpi
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
for (i in 1:N){
  ts[,i] <- (ts[,i] - mean(ts[,i])) / sd(ts[,i])
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
y_index <- seq(N)[colnames(ts) %in% y_names]

'%not_in%' <- Negate('%in%')
x_index <- which(seq(N) %not_in% y_index)
slow_index_not_y <- slow_index[slow_index %not_in% y_index] # w.r.t ts object

Y <- ts[,y_index]
X <- ts[,x_index]
X_slow <- ts[,slow_index_not_y]
print(colnames(Y))
T <- dim(Y)[1]
y_name <- "aRUCPI_slow"
i_name <- "RUCBIR=ECI"
y <- Y[,2]

results_cpi = c()
for (k in 0:5) {
  errors_FAVAR_cv_cpi_k <- my_tsCV(y = y, forecastfunction = fore_FAVAR, y_name=y_name,
                                   h = 12, X = X, Y = Y, X_slow = X_slow, K = k + (k==0), 
                                   i_name = i_name, 
                                   initial = 200, use_VAR = (k==0))
  res_FAVAR_cpi_k <- sqrt(colMeans(errors_FAVAR_cv_cpi_k^2, na.rm = TRUE))
  results_cpi = c(results_cpi, sum(res_FAVAR_cpi_k))
}

df_cpi = data.frame(k=seq(0,5), sum_rmse=results_cpi)

ggplot(data=df_cpi, aes(x=k, y=sum_rmse)) + geom_point(shape=15, size=3) + geom_line() + labs(x='K', y='Sum of RMSE') +
  ggtitle("CPI forecast error, Russia")

# unemp
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


results_unemp = c()
for (k in 0:5) {
  errors_FAVAR_cv_unemp_k <- my_tsCV(y = y, forecastfunction = fore_FAVAR, y_name=y_name,
                                   h = 12, X = X, Y = Y, X_slow = X_slow, K = k + (k==0), 
                                   i_name = i_name, 
                                   initial = 200, use_VAR = (k==0))
  res_FAVAR_unemp_k <- sqrt(colMeans(errors_FAVAR_cv_unemp_k^2, na.rm = TRUE))
  results_unemp = c(results_unemp, sum(res_FAVAR_unemp_k))
}

df_unemp = data.frame(k=seq(0,5), sum_rmse=results_unemp)

ggplot(data=df_unemp, aes(x=k, y=sum_rmse)) + geom_point(shape=15, size=3) + geom_line() + labs(x='K', y='Sum of RMSE') +
  ggtitle("Unemployment rate forecast error, Russia")


