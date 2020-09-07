function []=mymonthlybars(y,stdy,maxvalue,delta)
    
%********************************************
%Programa MYMONTHLYBARS.m: Este programa me obtiene un plot de monthly
%barras para 1 variable.
%
%Use: mymonthlybars(y,stdy,maxvalue,delta)
%
%donde:
% y: Vector y de 12 valores (monthly means).
% stdy: Vector std(y) de 12 valores (monthly means).
% maxvalue: Limite superior que deseo para el grafico.
% delta: distancia entre las lineas de la "grid on" que deseo.
%********************************************
x=[1:length(y)];
z1 = y + stdy;
z2 = y - stdy;

a1=area(x,z1);
grid on
set(gca,'Layer','top')
hold on
a2=area(x,z2);
grid on
set(gca,'Layer','top')
hold on
hb=bar(y(:),0.5,'grouped');
hold off
shading faceted
%color=[0 0 0.8]
color=[0 0 0]
set(hb(1),'Facecolor',color)
meses=['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'Xlim',[0 13],'XTick',[1:1:12],'XTickLabel',meses,'Ylim',[0 maxvalue],'YTick',[0:delta:maxvalue])
%..............................................
% $$$ set(gca,'Ycolor',[0 0 0],'FontSize',[5]);
%..............................................
set(gca,'Ycolor',[0 0 0],'FontSize',[12]);
%..............................................
set(a1,'LineStyle','none','FaceColor',[0.9 0.9 0.9]);
set(a2,'LineStyle','none','FaceColor',[1 1 1]);
grid on
