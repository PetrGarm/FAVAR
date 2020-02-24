library("tidyverse")
library("forecast")
library("vars")
library("lmtest")

df <- read.table("C:/Users/petrg/Desktop/ЦМП/BBE_Ddisk/BBE_Ddisk/2step/nsbalpanel.txt")
ts <- ts(df, start = c(1959,1), freq = 12)

# Standardize for some reason (correct PCA extraction probably)
for (i in 1:120){
  ts[,i] <- (ts[,i] - mean(ts[,i]))/sd(ts[,i])
}

ggAcf(ts[,1])
ggPacf(ts[,1])
ts.plot(ts[,1])

# x_index <- c(78,81,92,96,97,98,74,102,17,49,32,46,54,62,66,119,120)
slow_index <- c(1:53, 103:119)
y_index <- c(16, 108, 77)

'%not_in%' <- Negate('%in%')
x_index <- which(seq(1,120) %not_in% y_index)

slow_index_not_y <- which(slow_index %not_in% y_index)


Y <- ts[,y_index]
X <- ts[,x_index]
X_slow <- ts[,slow_index_not_y]
colnames(Y) <- c("IP", "CPI", "FFR")
ts.plot(Y[,3])


# number of latent factors
K <- 3
# number of Y
M <- length(y_index)

# PCA
prin_comp <- princomp(ts)
PCs <- prin_comp$scores[,seq(1, K)]
Fs <- prin_comp$scores[,seq(1, K)]
colnames(Fs) <- c('fac1', 'fac2', 'fac3')
PCs_slow <- princomp(X_slow)$scores[,seq(1, K)]

new_df = cbind(Y, PCs_slow)

for(i in 1:K){
  temp = cbind(PCs[,i], new_df)
  colnames(temp)[1] = "PC_old"
  lr <- lm(PC_old ~ ., data = temp)
  
  C_hat <- lr$fitted.values
  Y_coefs <- lr$coefficients[c(2, 3 ,4)]
  
  F_hat <- C_hat - Y_coefs[1] * Y[,1] - Y_coefs[2] * Y[,2] - Y_coefs[3] * Y[,3] 
  Fs[,i] <- F_hat
}

Y_for_VAR <- cbind(Fs,Y)
VARselect(Y_for_VAR)
var <- VAR(Y_for_VAR, p = 10)
summary(var)

causality(var, 'Fs.fac1')
causality(var, 'Fs.fac2')
causality(var, 'Fs.fac3')

y_hat = forecast(var, h=12)

feir <- irf(var, impulse = "Y.FFR", response = "Y.IP",
            n.ahead = 8, ortho = FALSE, runs = 1000)

plot(feir)

oir <- irf(var, impulse = "Y.IP", response = "Y.FFR",
            n.ahead = 8, ortho = TRUE, runs = 1000, seed = 322)

plot(oir)

