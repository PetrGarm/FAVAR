library("tidyverse")
library("forecast")
library("vars")
library("lmtest")

fore_FAVAR <- function(X, Y, K, y_name, h, y, use_VAR = F){
  start <- start(X)
  frequency <- frequency(X)
  X = X[1:length(y),]
  Y = Y[1:length(y),,drop=FALSE]
  
  N = dim(X)[2]
  M = dim(Y)[2]
  T = dim(X)[1]
  
  # extract PC from X
  XtX = t(X) %*% X
  spectral_decomposition <- eigen(XtX)
  eig_vec <- spectral_decomposition$vectors
  eig_val <- spectral_decomposition$values
  
  lam <- 1/sqrt(N) * eig_vec[,(1:K)]
  Fr <- data.matrix(X) %*% lam
  
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