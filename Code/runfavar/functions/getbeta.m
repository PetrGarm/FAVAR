function out=getbeta(xs,xn)
out=[];
for i=1:cols(xs)
  ytemp=xn(:,i);
  xtemp=[xs(:,i) ones(rows(xs),1)];
  bn=xtemp\ytemp;
  out=[out;bn(1)];
end
  