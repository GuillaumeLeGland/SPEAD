function [myblue,myblueinv]=mycolorbarblue()

x=[1,56];
xi=[1:56];
fila1=[1,1,1]; %white
fila2=[0,0,1]; %blue
J1=[fila1(1);fila2(1)]; %col1
J2=[fila1(2);fila2(2)]; %col2
J3=[fila1(3);fila2(3)]; %col3

J1I=interp1(x,J1,xi)';
J2I=interp1(x,J2,xi)';
J3I=interp1(x,J3,xi)';
blue=[J1I,J2I,J3I]; %(48,3)

%ANADO LOS EXTREMOS ROJO-FUERTE Y AZUL-FUERTE (los cojo del "jet" colormap):
pal=[0.9375,0.8750,0.8125,0.7500,0.6875,0.6250,0.5625,0.5000]';
ceros=zeros(8,1);
MAPazul=[ceros,ceros,pal]; %(8 x 3)
myblue=[blue;MAPazul]; %(64 x 3) donde "negativo-blanco, positivo-azul".

%OBTENGO TAMBIEN LA PALETA DE COLORES INVERSA:
myblueinv=flipud(myblue); %(64 x 3) donde "negativo-azul, positivo-blanco".

%............
% $$$ X=[-10:10];
% $$$ figure(1),plot(X),colormap(myblue),colorbar('horiz')
% $$$ figure(2),plot(X),colormap(myblueinv),colorbar('horiz')
% $$$ pause
%............
