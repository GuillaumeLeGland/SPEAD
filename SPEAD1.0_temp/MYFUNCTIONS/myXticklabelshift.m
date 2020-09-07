function [hplot] = myXticklabelshift(x,y)

monthlim1 = [0:30:360];
monthlim2 = 0.5*diff(monthlim1)+monthlim1(1:end-1)
xticklabels = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];

ax(1) = newplot;
set(gcf,'nextplot','add');
set(ax(1),'Xlim',[0 360],'XTick',monthlim2,'XTicklabel',xticklabels);
ax(2) = axes('position',get(ax(1),'position'),'Visible', 'off');
hplot = plot(x,y);
set(ax(2),'Xlim',[0 360],'XTick',monthlim1,'XTicklabel','');
grid on

