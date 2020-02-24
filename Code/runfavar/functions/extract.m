function [fac,lam]=extract(data,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function extracts the first k principal components from a t*n    %
%matrix 'data' and returns the factors (fac, t*k) as well as the       % 
%normalised 'loading' matrix (lam, n*k) - i.e. the matrix composed of  %
%eigenvectors of the data's covariance matrix such that the principal  %
%components are equal to:  fac = X*lam                                 %
%                                                                      % 
%Note 1: We normalise to ensure that lam'lam/n=I 
%Note 2: We normalise to ensure that lam'lam/n=I 
%. The loading matrix is normalised, so that      %
%or matrix i.e. for given  x=(a_{1}; a_{2}; ...; a_{n}) (where a_{i}   %
%are row vectors) it returns the vector or matrix:                     %
%                                                                      % 
%xld=(log(a_{p+1})-log(a_{1});log(a_{p+2})-log(a_{2}); ...             % 
%                                     ... ;log(a_{n})-log(a_{n-p}))    %       
%                                                                      % 
%Note: The function does not check whether the means of the variables  % 
%      are zero, if they are not it will return incorrect results      % 
%                                                                      % 
%Version is 1.01                                                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute the size of the data matrix
[t,n]=size(data);
%Assuming the series in the columns have empirical mean = 0, the 
%empirical variance-covariance matrix is given by 
xx=data'*data;
%Compute the eigenvector/eigenvalue decomposition of X'X
[evec,eval]=eig(xx);

%Order eigenvalues in descending order and save the sorting permutation in index
[eval,index]=sort(diag(eval),'descend');

%Reorder the eigenvectors accordingly (i.e. from those corresponding to 
%the largest eigenvalue to those corresponding to the smallest eigenvalue);
evc=zeros(n,n);
for i=1:n
   evc(:,i)=evec(:,index(i));
end

%Take the first k eigenvectors and normalise by square root of #data series
lam = sqrt(n)*evc(:,1:k);

%Compute the principal component vectors from standard formula; Note here
%all the principal componens are 'shortened' by the division 
fac=data*lam/n;

%Not sure why we bother normalising, though.

%Revision History:
%v 1.01
%Replaced:
%function [fac,lam]=extract1(data,k)
%by
%function [fac,lam]=extract(data,k)
%%
%Replace
%[eval,index]=sort(diag(eval));
%index=flipud(index); 		   %to get descending order
%by
%[eval,index]=sort(diag(eval),'descend');
