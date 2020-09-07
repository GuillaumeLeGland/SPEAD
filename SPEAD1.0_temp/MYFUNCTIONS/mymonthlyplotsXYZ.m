function []=mymonthlyplotsXYZ(x,y,z,xstd,ystd,zstd,maxvalue,delta)

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
%********************************************

%...............................
% $$$ maxvalue=max([x(:);y(:);z(:)]);
% $$$ delta=maxvalue/10;
%...............................
hp=plot([1:12],x(:),'k-',[1:12],y(:),'k--',[1:12],z(:),'k:');
hold on
he1=errorbar([1:12],x,xstd);
hold on
he2=errorbar([1:12],y,ystd);
hold on
he3=errorbar([1:12],z,zstd);
hold off
meses=['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
set(gca,'Xlim',[0 13],'XTick',[1:1:12],'XTickLabel',meses,'Ylim',[0 maxvalue],'YTick',[0:delta:maxvalue]);
set(hp(1),'Marker','.','MarkerEdgeColor','k','Markersize',10,'linestyle','-');
set(hp(2),'Marker','o','MarkerEdgeColor','k','Markersize',4,'linestyle','--');   
set(hp(3),'Marker','x','MarkerEdgeColor','k','Markersize',8,'linestyle',':'); 
%...
set(he1(1),'Color',[0 0 0],'linestyle','-');
set(he2(1),'Color',[0 0 0],'linestyle','--');
set(he3(1),'Color',[0 0 0],'linestyle',':');
%...
% $$$ set(he1(1),'Color',[0 0 0],'linestyle','-');
% $$$ set(he2(1),'Color',[0 0 0],'linestyle','-');
% $$$ set(he3(1),'Color',[0 0 0],'linestyle','-');
%...
set(he1(2),'Marker','.','MarkerEdgeColor','k','Markersize',10,'linestyle','none');
set(he2(2),'Marker','o','MarkerEdgeColor','k','Markersize',4,'linestyle','none');
set(he3(2),'Marker','x','MarkerEdgeColor','k','Markersize',4,'linestyle','none');
%..............................................
% $$$ set(gca,'Ycolor',[0 0 0],'FontSize',[6]); 
%..............................................
set(gca,'Ycolor',[0 0 0],'FontSize',[12]); 
%..............................................
%grid on
