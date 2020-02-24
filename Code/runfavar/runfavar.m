clear
addpath('functions');
[ data0 junk ]=xlsread('\data\datain.xls');
[ junk names ]=xlsread('\data\names.xls');

index=xlsread('\data\index.xls');
dindex=index(:,1); %dindex=1 for series that are log differenced dindex=3 differencing without logs
index=index(:,2);  %index=1 for 'fast moving' series

reps=5000; %replications
burn=4000; %burn in 
horizon=40;

%first difference the data where appropriate
data00=diffx(data0,dindex);
%standardise the data
data=standardise(data00);
betamat=getbeta(data,data00);%scaling factors to convert IRFs into original data units


%load policy rate and standardize it
z=xlsread('\data\baserate.xls');
z=z(2:end);
z=standardise(z);



KK=2;  %number of factors
L=4;  %number of lags in the VAR
N=KK+1; %number of Variables in var K factors plus the interest rate
NN=cols(data);% size of the panel
T=rows(data);
%step 1 of the algorithm set starting values and priors

%get an intial guess for the factor via principal components
pmat=extract(data,KK);
beta0=[pmat(1,:) z(1)];  %state vector S[t-1/t-1]
for i=1:L-1
    beta0=[beta0 zeros(1,N)];
end
ns=cols(beta0);
P00=eye(ns);  %P[t-1/t-1]
rmatin=ones(NN,1); %arbitrary starting value for the variance of the idiosyncratic component
Sigmain=eye(N);  %arbitrary starting value for the variance of VAR errors

%prior for the factor loadings,variances and VAR
fload0=zeros(KK+1,1);
vfload0=eye(KK+1);
vg=0.01;
tg=5;
LAMDAP=0.1; %tightness of prior smaller=tighter
EPSILON=1/1000; %tightness of prior on constant, small=looser
[yd,xd]=getdummies([pmat z],LAMDAP,EPSILON,L);

mm=1;
for m=1:reps;

%gibbs sampling 
 disp(strcat('REPS=',num2str(m)));

%step 2 sample factor loadings
[fload,floadr,rmat]=getobseqparameters(data,pmat,z,fload0,vfload0,vg,tg,rmatin,index);
rmatin=rmat;

Y=[pmat z];
X=[];
for i=1:L
    X=[X lag0(Y,i)];
end
X=[X ones(rows(Y),1)];
Y=Y(L+1:end,:);
X=X(L+1:end,:);
Y0=[Y;yd];
X0=[X;xd];
 [beta,Sigma]=drawvar(Y0,X0,Sigmain,N,L);
 Sigmain=Sigma;

beta1=reshape(beta,N*L+1,N);



%step 5 prepare matrices for the state space
%Y=H*factors+e
%factors=MU+F*Factors(-1)+v
%e~N(0,R)
%v~N(0,Q)

%matrix of factor loadings
H=zeros(NN,N*L);
H(1:rows(fload),1:KK+1)=[fload floadr];
H(rows(floadr)+1,KK+1)=1;
%matrix R
R=diag([rmat;0]);
%vector MU
MU=[beta1(end,:)';zeros(N*(L-1),1)]';
%matrix F
F=[beta1(1:N*L,:)';eye(N*(L-1),N*L)];
%matrix Q
Q=zeros(rows(F),rows(F));
Q(1:N,1:N)=Sigma;


%Note, Lines 105 and 107 call the Carter Kohn code written by Bernanke
%Et.al (2005). To use this, comment in and comment out 110 to 167
% indexnM=find(repmat(1:KK+1,L,1)'<KK+1); %index of factors and their lags
% [Sdraw,S]=kfgibbsnv([data z],beta0',P00,H(:,1:KK+1),R,F,Q,cols(z),indexnM); %use bernanke et.al's code
% pmat=Sdraw(:,1:KK);


%Carter and Kohn algorithm to draw the factor
beta_tt=zeros(T,ns);          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
% %%%%%%%%%%%Step 6a run Kalman Filter
i=1;
x=H;
%Prediction
beta10=MU+beta0*F';
p10=F*P00*F'+Q;
yhat=(x*(beta10)')';                                                
eta=[data(i,:) z(i,:)]-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*invpd(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
beta_tt(i,:)=beta11;
ptt(i,:,:)=p11;
for i=2:T
    %Prediction
beta10=MU+beta11*F';
p10=F*p11*F'+Q;
yhat=(x*(beta10)')';                                               
eta=[data(i,:) z(i,:)]-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*invpd(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
ptt(i,:,:)=p11;
beta_tt(i,:)=beta11;
end
% Backward recursion to calculate the mean and variance of the distribution of the state
%vector
beta2 = zeros(T,ns);   %this will hold the draw of the state variable
jv1=1:KK; %index of state variables to extract
jv=1:N;
wa=randn(T,ns);
f=F(jv,:);
q=Q(jv,jv);
mu=MU(jv);
i=T;  %period t
p00=squeeze(ptt(i,jv1,jv1)); 
beta2(i,:)=beta_tt(i,:);
beta2(i,jv1)=beta_tt(i:i,jv1)+(wa(i:i,jv1)*cholx(p00));   %draw for beta in period t from N(beta_tt,ptt)
%periods t-1..to .1
for i=T-1:-1:1
pt=squeeze(ptt(i,:,:));
bm=beta_tt(i:i,:)+(pt*f'*inv(f*pt*f'+q)*(beta2(i+1:i+1,jv)-mu-beta_tt(i,:)*f')')';  
pm=pt-pt*f'*inv(f*pt*f'+q)*f*pt;  
beta2(i,:)=bm;

beta2(i:i,jv1)=bm(jv1)+(wa(i:i,jv1)*cholx(pm(jv1,jv1)));  

end


pmat=beta2(:,1:KK);   %update the factors

if m>burn
    %compute impulse response
    A0=cholx(Sigma);
    shock=zeros(1,N);
    shock(end)=1;
    yhat=irfsim(beta1,N,L,A0,shock,horizon);
   
yhat1=yhat*H(:,1:KK+1)';  %impulse response for the panel


irfmat(mm,1:horizon,1:NN+1)=(yhat1);

mm=mm+1;
end
    

end



irf=prctile(irfmat,[50 16 84],1);


figure(1)
j=1
for i=1:size(irf,3)-1
subplot(4,10,j)
plotx1(squeeze(irf(:,:,i))'.*betamat(i));
title(strcat('\fontsize{8}', names(i)))
j=j+1
axis tight
end