function plotx1(y)
set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;1 0 0;0 0 1]);
cu=y(:,2);
cl=y(:,3);

t=size(y,1);
h=0:t-1;
h=h';


% hhx=fill([h(2); h(1:end); flipud([h(1:end); h(end)])],[cu1(1); cl1(1:end); flipud([cu1(1:end); cl1(size(cl1,1))])],'b');
% set(hhx,'edgecolor',[1 1 1]);
% set(hhx,'facecolor',[0.90 0.90 1]);
% 
% hold on;

hh=fill([h(2); h(1:end); flipud([h(1:end); h(end)])],[cu(1); cl(1:end); flipud([cu(1:end); cl(size(cl,1))])],'b');
set(hh,'edgecolor',[0.75 0.75 1]);
set(hh,'facecolor',[0.75 0.75 1]);

hold on;

plot(h,y(:,1),'LineWidth',2);

% hold on;
% 
% cu=y(:,5);
% cl=y(:,6);
% 
% hh=fill([h(2); h(1:end); flipud([h(1:end); h(end)])],[cu(1); cl(1:end); flipud([cu(1:end); cl(size(cl,1))])],'b');
% set(hh,'edgecolor',[0.95 0.75 1]);
% set(hh,'facecolor',[0.95 0.75 1]);
% 
% hold on;
% 
% plot(h,y(:,4),'LineWidth',2);
% 

