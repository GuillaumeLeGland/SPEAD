function []=mymonthlybarsXY(x,y,maxvalue,delta,type)
    
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
% maxvalue: Limite superior que deseo para el grafico.
% delta: distancia entre las lineas de la "grid on" que deseo.
% type: 'SIcorrcoef' o 'NOcorrcoef' (segun si quiero o no calcular el coef.correlation entre x e y)
%********************************************
%...............................
% $$$ maxvalue=max([x(:);y(:)]);
% $$$ delta=maxvalue/10;
%...............................
hb=bar([x(:),y(:)],1.0,'grouped');
shading faceted
set(hb(1),'Facecolor',[0 0 0.8])
set(hb(2),'Facecolor',[0 0.8 0])
meses=['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'Xlim',[0 13],'XTick',[1:1:12],'XTickLabel',meses,'Ylim',[0 maxvalue],'YTick',[0:delta:maxvalue]);
%..............................................
% $$$ set(gca,'Ycolor',[0 0 0],'FontSize',[6]); 
%..............................................
set(gca,'Ycolor',[0 0 0],'FontSize',[10]); 
%..............................................
grid on
if strcmp(type,'SIcorrcoef')
    rho=mycorrcoef(x(:),y(:),'Spearman');
    corr=num2str(round(100*rho)/100);
    htxt1=text(0.5,0.93*maxvalue,['\rho = ',corr]);
    %color=[0 0 0.8]; %azul.
    color=[0.9 0 0]; %rojo
    %.......................................
% $$$     set(htxt1,'Fontsize',[6],'Color',color)
    %.......................................
    set(htxt1,'Fontsize',[8],'Color',color)
    %.......................................
else strcmp(type,'NOcorrcoef')
end
