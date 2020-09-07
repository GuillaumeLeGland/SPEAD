%function []=mymonthlyXticklabels(y,A)
function []=mymonthlyXticklabels()
close all
clear all

monthlim=[1,31,59,90,120,151,181,212,243,273,304,334,365];
monthlim2=0.5*diff(monthlim)+monthlim(1:end-1);
xticklabels=['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
 
y=rand(1,365)*100;
A=rand(150,365)*100;

figure(1)
hp1=plot([1:365],y,'-r.');
axis square
set(gca,'Xlim',[1,365],'XTick',monthlim2,'XTickLabel',xticklabels);
set(gca,'Ylim',[0,100],'YTick',[0:10:100],'YTickLabel',[0:10:100]);
original = gca
position = get(original, 'Position')
temporary = axes('Position', position,'Visible', 'off')
hp1=plot([1:365],y,'-b.');
axis square
set(temporary,'Xlim',[1,365],'XTick',monthlim,'XTicklabel','','Yticklabel','');
grid on

figure(1)
hp1=plot([1:365],y,'-r.');
hleg1=legend(hp1,'Legenda A',1)
axis square
set(gca,'Xlim',[1,365],'XTick',monthlim2,'XTickLabel',xticklabels);
set(gca,'Ylim',[0,100],'YTick',[0:10:100],'YTickLabel',[0:10:100]);
original = gca
position = get(original, 'Position')
temporary = axes('Position', position,'Visible', 'off')
hp1=plot([1:365],y,'-b.');
hleg2=legend(hp1,'Legenda B',2)
axis square
set(temporary,'Xlim',[1,365],'XTick',monthlim,'XTicklabel','','Yticklabel','');
grid on


figure(2)
himg=imagesc(A)
hc=colorbar('vertic');
colormap('cool')
set(gca,'Xlim',[1,365],'XTick',monthlim2,'XTickLabel',xticklabels)
set(gca,'Ylim',[0,150],'YTick',[0:25:150],'Yticklabel',[0:25:150]);
set(gca,'Fontsize',[4])
set(hc,'Fontsize',[4])
original = gca
position = get(original, 'Position')
temporary = axes('Position', position,'Visible', 'off')
himg=imagesc(A)
set(temporary,'Xlim',[1,365],'XTick',monthlim,'XTicklabel','');
set(temporary,'Ylim',[0,150],'YTick',[0:25:150],'Yticklabel','');
grid on
pause

%******************************************************
% $$$ y=rand(1,365);
% $$$ figure(1)
% $$$ ax(1)=newplot;
% $$$ set(gcf,'nextplot','add');
% $$$ set(ax(1),'Xlim',[1,365],'XTick',monthlim2,'XTicklabel',xticklabels);
% $$$ ax(2)=axes('position',get(ax(1),'position'),'Visible', 'off');
% $$$ plot(1:365,y)
% $$$ set(ax(2),'Xlim',[1,365],'XTick',monthlim,'XTicklabel','');
% $$$ grid on
% $$$ 
% $$$ y=rand(1,365)*100;
% $$$ figure(10)
% $$$ ax(1)=newplot;
% $$$ set(gcf,'nextplot','add');
% $$$ set(ax(1),'Xlim',[1,365],'XTick',monthlim2,'XTicklabel',xticklabels,'Ylim',[0,100],'YTick',[0:10:100]);
% $$$ ax(2)=axes('position',get(ax(1),'position'),'Visible', 'off');
% $$$ plot(1:365,y)
% $$$ set(ax(2),'Xlim',[1,365],'XTick',monthlim,'XTicklabel','','Yticklabel','');
% $$$ grid on
% $$$ 
% $$$ y=rand(1,365)*100;
% $$$ figure(100)
% $$$ hp1=plot([1:365],y,'-r.');
% $$$ set(gca,'Xlim',[1,365],'XTick',monthlim2,'XTickLabel',xticklabels,'Ylim',[0,100],'YTick',[0:10:100]);
% $$$ original = gca
% $$$ position = get(original, 'Position')
% $$$ temporary = axes('Position', position,'Visible', 'off')
% $$$ hp1=plot([1:365],y,'-r.');
% $$$ set(temporary,'Xlim',[1,365],'XTick',monthlim,'XTicklabel','','Yticklabel','');
% $$$ grid on

%A=rand(150,365)*10;
figure(2)
subplot(2,2,1)
himg=imagesc(A)
hc=colorbar('horiz');
set(gca,'Xlim',[1,365],'XTick',monthlim2,'XTickLabel',xticklabels)
set(gca,'Ylim',[0,150],'YTick',[0:25:150],'Yticklabel',[0:25:150]);
set(gca,'Fontsize',[4])
set(hc,'Fontsize',[4])
original = gca
position = get(original, 'Position')
temporary = axes('Position', position,'Visible', 'off')
himg=imagesc(A)
set(temporary,'Xlim',[1,365],'XTick',monthlim,'XTicklabel','');
set(temporary,'Ylim',[0,150],'YTick',[0:25:150],'Yticklabel','');
grid on

subplot(2,2,2)
% $$$ himg=imagesc(A)
set(gca,'Xlim',[1,365],'XTick',monthlim2,'XTickLabel',xticklabels)
set(gca,'Ylim',[0,150],'YTick',[0:25:150],'Yticklabel',[0:25:150]);
set(gca,'Fontsize',[4])
set(hc,'Fontsize',[4])
original = gca
position = get(original, 'Position')
temporary = axes('Position', position,'Visible', 'off')
himg=imagesc(A)
hc=colorbar('horiz');
set(temporary,'Xlim',[1,365],'XTick',monthlim,'XTicklabel','');
set(temporary,'Ylim',[0,150],'YTick',[0:25:150],'Yticklabel','');
grid on

subplot(2,2,3)
himg=imagesc(A)
colorbar('horiz')
set(gca,'Xlim',[1,365],'XTick',monthlim2,'XTickLabel',xticklabels,'Ylim',[0,150],'YTick',[0:25:150]);
original = gca
position = get(original, 'Position')
temporary = axes('Position', position,'Visible', 'off')
% $$$ himg=imagesc(A)
% $$$ colorbar('horiz')
set(temporary,'Xlim',[1,365],'XTick',monthlim,'XTicklabel','','Yticklabel','');
grid on

%********************************
% $$$ figure(2)
% $$$ ax(1)=newplot;
% $$$ set(gcf,'nextplot','add');
% $$$ set(ax(1),'Xlim',[1,365],'XTick',monthlim,'XTicklabel',monthlim);
% $$$ ax(2)=axes('position',get(ax(1),'position'));
% $$$ plot(1:365,y)
% $$$ set(ax(2),'XAxisLocation','top');
% $$$ set(ax(2),'Xlim',[1,365],'XTick',monthlim2,'XTicklabel',xticklabels);

