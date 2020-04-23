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

slow_index_not_y <- slow_index[slow_index %not_in% y_index] # w.r.t ts object


Y <- ts[,y_index]
X <- ts[,x_index]
X_slow <- ts[,slow_index_not_y]
colnames(Y) <- c("IP", "CPI", "FFR")
ts.plot(Y[,3])


K = 3 # Number of factors
lags = 13 # in VAR
N = dim(X)[2]
M = dim(Y)[2]
T = dim(X)[1]

# extract PC from X
XtX = t(X) %*% X
spectral_decomposition <- eigen(XtX)
eig_vec <- spectral_decomposition$vectors
eig_val <- spectral_decomposition$values

lam <- 1/sqrt(N) * eig_vec[,(1:K)]
F0 <- data.matrix(X) %*% lam 

X_slow <- ts[,slow_index_not_y]

XtX = t(X_slow) %*% X_slow
spectral_decomposition <- eigen(XtX)
eig_vec <- spectral_decomposition$vectors
eig_val <- spectral_decomposition$values

nslow <- dim(X_slow)[2]
lam <- 1 / sqrt(nslow) * eig_vec[,(1:K)]
Fslow0 <- data.matrix(X_slow) %*% lam 

# Factors rotation 
ones_vec = rep(1, T)
Ffast = data.matrix(Y[,3])
k1 = 1 # ? ?? ?
ly = cbind(ones_vec, Ffast, Fslow0)

svd <- svd(ly)
vl <-svd$u
d <- 1 / svd$d
vr <- svd$v

b = (vr * matrix(data = rep(d, dim(vr)[1]), ncol = dim(vr)[1], byrow = TRUE)) %*% (t(vl) %*% F0)

Ffast_ktimes = matrix(rep(Ffast, K), ncol = K)
coeffs = matrix(rep(coeffs, dim(Ffast_ktimes)[1]), ncol=K, byrow = TRUE)
Fr = F0 - Ffast_ktimes*coeffs

# FAVAR with clean factors
Y_for_VAR = cbind(Fr, Y)
VARselect(Y_for_VAR)
p = 3 

var <- VAR(Y_for_VAR, p)
summary(var)

causality(var, 'Fr.Series.1')
causality(var, 'Fr.Series.2')
causality(var, 'Fr.Series.3')

y_hat = forecast(var, h=12)
autoplot(y_hat$forecast$Y.IP,) + ggtitle("US IP forecast from FAVAR(3)[k=3]") + ylab("US IP")
