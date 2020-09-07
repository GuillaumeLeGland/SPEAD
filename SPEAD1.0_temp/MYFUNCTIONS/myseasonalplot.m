function []=myseasonalplot(y,ymin,ymax,dy)

monthlim=[1,31,59,90,120,151,181,212,243,273,304,334,365];
monthlim2=0.5*diff(monthlim)+monthlim(1:end-1);
xticklabels=['J','F','M','A','M','J','J','A','S','O','N','D','J']';

%%%%%%%%%%
%1st AXES:
%%%%%%%%%%
hp=plot([1:365],y,'-r.');
%...............
%Ha de estar before de "colorbar".
ax1 = gca;
ax1position = get(ax1, 'Position');
%...............
set(ax1,'Xlim',[1 365],'XTick',monthlim2,'XTickLabel',xticklabels)
set(ax1,'Ylim',[ymin ymax],'YTick',[ymin:dy:ymax],'YTickLabel',[ymin:dy:ymax]);
%axis square 

%%%%%%%%%%
%2nd AXES:
%%%%%%%%%%
ax2 = axes('Position', ax1position,'Visible', 'off');
hp=plot([1:365],y,'-r.');
set(ax2,'Xlim',[1 365],'XTick',monthlim,'XTickLabel','')
set(ax2,'Ylim',[ymin ymax],'YTick',[ymin:dy:ymax],'YTickLabel',[ymin:dy:ymax]);
set(gca,'Fontsize',10)
%axis square 
grid on
