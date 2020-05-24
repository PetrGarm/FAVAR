library("tidyverse")
library("forecast")
library("vars")
library("lmtest")

fore_FAVAR <- function(X, Y, X_slow, K, y_name, i_name, h, y, use_VAR = F){
  start <- start(X)
  frequency <- frequency(X)
  X = X[1:length(y),]
  Y = Y[1:length(y),]
  X_slow = X_slow[1:length(y),]
  
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
  
  #extract PC from X_slow
  XtX = t(X_slow) %*% X_slow
  spectral_decomposition <- eigen(XtX)
  eig_vec <- spectral_decomposition$vectors
  eig_val <- spectral_decomposition$values
  
  nslow <- dim(X_slow)[2]
  lam <- 1 / sqrt(nslow) * eig_vec[,(1:K)]
  Fslow0 <- data.matrix(X_slow) %*% lam 
  
  # Factors rotation 
  ones_vec = rep(1, T)
  Ffast = data.matrix(Y[,i_name]) # interest rate
  k1 = dim(Ffast)[2] # (num of Ffast)
  ly = cbind(ones_vec, Ffast, Fslow0)
  
  svd <- svd(ly)
  vl <-svd$u
  d <- 1 / svd$d
  vr <- svd$v
  
  b = (vr * matrix(data = rep(d, dim(vr)[1]), ncol = dim(vr)[1], byrow = TRUE)) %*% (t(vl) %*% F0)
  
  Ffast_ktimes = matrix(rep(Ffast, K), ncol = K)
  coeffs = matrix(rep(b[2:(k1+1),], dim(Ffast_ktimes)[1]), ncol=K, byrow = TRUE)
  Fr = F0 - Ffast_ktimes*coeffs
  
  # FAVAR with clean factors
  if(use_VAR == F){
    Y_for_VAR = cbind(Fr, Y)
    factor_names = paste('factor', 1:K)
    colnames(Y_for_VAR)[1:K] <- factor_names
  }
  else 
    Y_for_VAR = Y
  
  Y_for_VAR <- ts(Y_for_VAR, start = start, frequency = frequency)
  VARselect(Y_for_VAR)
  p = VARselect(Y_for_VAR)$selection[1] # according AIC 
  
  var <- VAR(Y_for_VAR, p)
  y_hat <- forecast::forecast(var, h=h)
  return(y_hat)
}