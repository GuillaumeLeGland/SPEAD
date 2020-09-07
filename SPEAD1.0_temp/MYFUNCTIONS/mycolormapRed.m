function [myred,myredinv]=mycolormapRed()

x=[1,56];
xi=[1:56];
fila1=[1,1,1]; %white
fila2=[1,0,0]; %red
J1=[fila1(1);fila2(1)]; %col1
J2=[fila1(2);fila2(2)]; %col2
J3=[fila1(3);fila2(3)]; %col3

J1I=interp1(x,J1,xi)';
J2I=interp1(x,J2,xi)';
J3I=interp1(x,J3,xi)';
red=[J1I,J2I,J3I]; %(48,3)

%ANADO EL EXTREMO RED-FUERTE (los cojo del "jet" colormap):
pal=[0.9375,0.8750,0.8125,0.7500,0.6875,0.6250,0.5625,0.5000]';
ceros=zeros(8,1);
MAProjo=[pal,ceros,ceros]; %(8 x 3)
myred=[red;MAProjo]; %(64 x 3) donde "negativo-blanco, positivo-red".

%OBTENGO TAMBIEN LA PALETA DE COLORES INVERSA:
myredinv=flipud(myred); %(64 x 3) donde "negativo-red, positivo-blanco".

%............
% $$$ X=[-10:10];
% $$$ figure(1),plot(X),colormap(myred),colorbar('horiz')
% $$$ figure(2),plot(X),colormap(myredinv),colorbar('horiz')
% $$$ pause
%............
