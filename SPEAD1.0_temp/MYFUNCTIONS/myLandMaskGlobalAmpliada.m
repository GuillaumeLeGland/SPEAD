function [GLOBE,Landampliada]=myLandMaskGlobalAmpliada(Land,varargin)
%******************************************
%Use: [Landampliada]=myLandMaskGlobalAmpliada([winx,winy]);
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
if numel(varargin)==0
    winx = 3; %Lat: si h=3 se aumenta en un "grado" (ie. pixel 1x1) la linea de costa.
    winy = 3; %Lon: si h=3 se aumenta en un "grado" (ie. pixel 1x1) la linea de costa.
elseif numel(varargin)==1
    runwin = varargin{1};
    winx = runwin(1); %Lat.
    winy = runwin(2); %Lon.
end
conwin = ones(winx,winy);
LAND = zeros(180,360);
LAND(Land) = -1;  %pongo -1 en la "Land".
CL = conv2(LAND,conwin,'same');
Landampliada = find(CL < 0);

TERRA = ones(180,360)*(-1);
TERRA(Landampliada) = 0;
TERRA(Land) = +1;
Coastline = find(TERRA==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ME QUEDO SOLO CON LA LINEA DE COSTA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GLOBE = zeros(180,360);
GLOBE(Coastline) = 1;
cmap = jet;
figure(1)
subplot(2,2,1)
imagesc(TERRA)
mycolorbar('horiz')
colormap(cmap)
subplot(2,2,2)
imagesc(GLOBE)
mycolorbar('horiz')
colormap(cmap)
pause(1)
