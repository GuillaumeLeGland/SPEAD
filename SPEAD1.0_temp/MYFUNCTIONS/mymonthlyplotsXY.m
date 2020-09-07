function []=mymonthlyplotsXY(x,y,xstd,ystd,maxvalue,delta,type)

%********************************************
%Programa MYMONTHLYBARS.m: Este programa me obtiene un plot de monthly
%barras para 2 variables, y tb da el coeficiente de correlacion entre
%ellas.
%
%Use: mymonthlyplots(x,y,maxvalue,delta)
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
hp=plot([1:12],x(:),'k-',[1:12],y(:),'k--');
hold on
he1=errorbar([1:12],x,xstd);
hold on
he2=errorbar([1:12],y,ystd);
hold off
meses=['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'Xlim',[0 13],'XTick',[1:1:12],'XTickLabel',meses,'Ylim',[0 maxvalue],'YTick',[0:delta:maxvalue]);
set(hp(1),'Marker','.','MarkerEdgeColor','k','Markersize',10,'linestyle','-');
set(hp(2),'Marker','o','MarkerEdgeColor','k','Markersize',4,'linestyle','--');   
set(he1(1),'Color',[0 0 0],'linestyle','-');
set(he2(1),'Color',[0 0 0],'linestyle',':');
set(he1(2),'Marker','.','MarkerEdgeColor','k','Markersize',10,'linestyle','none');
set(he2(2),'Marker','o','MarkerEdgeColor','k','Markersize',4,'linestyle','none');
%..............................................
% $$$ set(gca,'Ycolor',[0 0 0],'FontSize',[6]); 
%..............................................
set(gca,'Ycolor',[0 0 0],'FontSize',[12]); 
%..............................................
%grid on
if strcmp(type,'SIcorrcoef')
    rho=mycorrcoef(x(:),y(:),'Spearman');
    corr=num2str(round(100*rho)/100);
    htxt1=text(0.5,0.93*maxvalue,['\rho = ',corr]);
    %color=[0 0 0.8]; %azul.
    %color=[0.9 0 0]; %rojo
    color=[0 0 0]; %negro.
    %.......................................
% $$$     set(htxt1,'Fontsize',[6],'Color',color)
    %.......................................
    set(htxt1,'Fontsize',[12],'Color',color)
    %.......................................
else strcmp(type,'NOcorrcoef')
end
