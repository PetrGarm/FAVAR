function y=standardise(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function inputs a matrix x and returns its 'standardised'        %
%equivalent - i.e. a matrix obtained by subtracting the mean of every  %
%column and dividing all column elements by the respective columns     %
% standard deviation.                                                  %
%                                                                      % 
%Version is 1.00                                                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=rows(x); n=cols(x);

%Create matrices m and s (of dimensions identical to x) which contain
%the means and standard deviations of every column in x
m=repmat(mean(x),t,1); s=repmat(std(x),t,1);

%Subtract the mean from all column elements and divide them by the standard deviation
y=(x-m)./s;