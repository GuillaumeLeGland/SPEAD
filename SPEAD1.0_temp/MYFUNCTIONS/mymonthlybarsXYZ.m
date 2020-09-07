function []=mymonthlybarsXYZ(x,y,z,maxvalue,delta)
    
%********************************************
%Programa MYMONTHLYBARS.m: Este programa me obtiene un plot de monthly
%barras para 2 variables, y tb da el coeficiente de correlacion entre
%ellas.
%
%Use: mymonthlybars(x,y,maxvalue,delta)
%
%donde:
% x: Vector x de 12 valores (monthly means).
% y: Vector y de 12 valores (monthly means).
% z: Vector z de 12 valores (monthly means).
% maxvalue: Limite superior que deseo para el grafico.
% delta: distancia entre las lineas de la "grid on" que deseo.
% type: 'SIcorrcoef' o 'NOcorrcoef' (segun si quiero o no calcular el coef.correlation entre x e y)
%********************************************
%...............................
% $$$ maxvalue=max([x(:);y(:);z(:)]);
% $$$ delta=maxvalue/10;
%...............................
hb=bar([x(:),y(:),z(:)],1.0,'grouped');
shading faceted
%................................
% $$$ set(hb(1),'Facecolor',[0 0 0.8])
% $$$ set(hb(2),'Facecolor',[0 0.8 0])
% $$$ set(hb(3),'Facecolor',[0.8 0 0])
%................................
set(hb(1),'Facecolor',[0.4 0.4 0.4])
set(hb(2),'Facecolor',[0.8 0.8 0.8])
set(hb(3),'Facecolor',[0 0 0])
%................................
meses=['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'Xlim',[0 13],'XTick',[1:1:12],'XTickLabel',meses,'Ylim',[0 maxvalue],'YTick',[0:delta:maxvalue]);
% $$$ set(gca,'Ycolor',[0 0 0],'FontSize',[4]); 
set(gca,'Ycolor',[0 0 0],'FontSize',[12]); 
grid on
