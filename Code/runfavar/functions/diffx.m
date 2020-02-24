function data=diffx(data0,dindex)
data=[];
for i=1:cols(data0);
    if dindex(i)==1
        dat=log(data0(:,i));
        dat=diff(dat)*100;
    elseif dindex(i)==3
        dat=diff(data0(:,i));
    else
        dat=data0(2:end,i);
    end
    data=[data dat];
end
