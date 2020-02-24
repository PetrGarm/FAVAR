function B=getreg(Y,X,B0,Sigma0,sigma2)

V=invpd(invpd(Sigma0)+(1/sigma2)*(X'*X));
M=V*(invpd(Sigma0)*B0+(1/sigma2)*X'*Y); 

B=M+(randn(1,cols(X))*chol(V))';


