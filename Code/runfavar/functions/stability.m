function S=stability(beta,n,l)

%coef   (n*l+1)xn matrix with the coef from the VAR
%l      number of lags
%n      number of endog variables
%FF     matrix with all coef
%S      dummy var: if equal one->stability
coef=reshape(beta,n*l+1,n);
%coef
%coef
FF=zeros(n*l,n*l);
FF(n+1:n*l,1:n*(l-1))=eye(n*(l-1),n*(l-1));

%Note that the coef will be transposed (CHECK THIS IS CORRECT)

for i=1:l
    %FF(1:n,1+n*(i-1):n+n*(i-1))=coef(1+n*(i-1):n+n*(i-1),1:n)';
    FF(1:n,1+n*(i-1):n+n*(i-1))=coef(1+n*(i-1):n+n*(i-1),1:n);
    %Bc{i}=eye(2)
end
%FF
%FF=FF';
ee=max(abs(eig(FF)));
S=ee>=1;


