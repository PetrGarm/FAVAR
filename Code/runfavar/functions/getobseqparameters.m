function [fload,floadr,rmat]=getobseqparameters(data,pmat,z,fload0,vfload0,vg,tg,rmatin,index)
NN=cols(data);
KK=cols(pmat);
fload=zeros(NN,KK);
floadr=zeros(NN,1);
rmat=zeros(NN,1);
temp=eye(KK);
temp1=zeros(KK,1);
k=0; %counts 
for i=1:KK
   y=data(:,i);
   if index(i)==0
        ff=temp(i,:)';
          x=pmat;
    else
        ff=[temp(i,:) temp1(i,:)]';
         x=[pmat z];
         k=k+1;
        
   end
    
   %save
    if index(i)==0;
        fload(i,:)=ff';
      
    else
            fload(i,1:KK)=ff(1:KK)';
        floadr(i,:)=ff(end);
    end
    error= y-x*ff;
    rmati= IG(tg,vg,error);
    rmat(i)=rmati;
end
   
 
   
for i=KK+1:NN
    
    if index(i)==0
        y=data(:,i);
        x=pmat;
        f0=fload0(1:KK,:);
        vf0=vfload0(1:KK,1:KK);
    else
        if k<=KK
        y=data(:,i); 
        x=pmat;
        f0=fload0(1:KK,:);
        vf0=vfload0(1:KK,1:KK);
        else
        y=data(:,i); 
        x=[pmat z];
        f0=fload0;
        vf0=vfload0;
        end
        k=k+1;
    end
   
    %draw
    ff=getreg(y,x,f0,vf0,rmatin(i));
    
    %save
    if index(i)==0;
        fload(i,:)=ff';
    else
        if k<=KK
            fload(i,1:KK)=ff(1:KK)';
        else
       fload(i,1:KK)=ff(1:KK)';
        floadr(i,:)=ff(end);
        end
    end
    error= y-x*ff;
    rmati= IG(tg,vg,error);
    rmat(i)=rmati;
end




