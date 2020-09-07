function [Landampliada,Coastline]=myglobelandampliada()
%******************************************
%Use: [Landampliada,Coastline]=corrGLOBAL_landampliada;
%******************************************
%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
%OBTENGO UN MAPA DE LAND+NO-DATA AMPLIADO (AUMENTO LA LINEA DE COSTA
%UNA DISTANCIA h DE LA WINDOW):
%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
%-----------------------------------------------------------------------
%NOTA: El metodo de la convolucion es como si en los pixels con -1 
%situo una ventana que abarcara el entorno del pixel y le asignara -1 
%a ese entorno tambien. Esto se hace para evitar coger ptos que
%queden tan cerca de Tierra que la window (usada en el calculo de la
%tempevol) toque la Tierra en algun momento (ponemos -888 en el 
%entorno de la Tierra con una ventana del mismo tamaño que la window 
%usada en la tempevol)
%-----------------------------------------------------------------------
[GLOBE,Land]=myglobeland;
h=3; %si h=3 se quita un "grado" (un pixel 1x1) de la linea de costa.
window=ones(h,h);
LAND=zeros(180,360);
LAND(Land)=-1;  %pongo -1 en la "Land".
CL=conv2(LAND,window,'same');
Landampliada=find(CL<0);

TERRA=ones(180,360)*(-1);
TERRA(Landampliada)=0;
TERRA(Land)=+1;
%..................
figure(1)
imagesc(TERRA)
colormap(jet)
colorbar('horiz')
%.................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ME QUEDO SOLO CON LA LINEA DE COSTA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Coastline=find(TERRA==0);
W=zeros(180,360);
W(Coastline)=1;
figure(2)
imagesc(W)
colormap(jet)
colorbar('horiz')
