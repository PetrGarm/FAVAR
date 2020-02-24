function [beta,sigma]=drawvar(Y0,X0,sigmain,N,L)
%conditional mean of the VAR coefficients
  mstar=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx=invpd(xx);  %inv(X0'X0) to be used later in the Gibbs sampling algorithm
  
   vstar=kron(sigmain,ixx);
       chck=-1;                 %make sure VAR is stationary
while chck<0
       beta=mstar+(randn(1,N*(N*L+1))*chol(vstar))';
S=stability(beta,N,L);
if S==0
    chck=10;
end
end
       
       %draw covariance
       e=Y0-X0*reshape(beta,N*L+1,N);
    scale=e'*e;
    sigma=iwpq(rows(Y0),invpd(scale));