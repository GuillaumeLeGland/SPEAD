function [jetflag,jetflaginv]=myjetflagcolorbar(varargin)

%*********************************************************************
%Use: [jetflag,jetflaginv]=myjetflagcolorbar(RedMax,BlueMax)
%
%Red Max: Value between 0.5 and 1.0 (strongnes of Red in RGB palete)
%Blue Max: Value between 0.5 and 1.0 (strongnes of Blue in RGB palete)
%*********************************************************************
    
%........................
nvarargin=length(varargin);
if nvarargin==2
    RedMax=varargin{1};
    BlueMax=varargin{1};
else
    RedMax=1.0;
    BlueMax=1.0;
end
%........................
if RedMax==1 & BlueMax==1
    x=[1,24,48];
    xi=[1:48];
else
    x=[1,32,64];
    xi=[1:64];
end
%........................
fila1=[RedMax,0,0]; %red
fila2=[1,1,1]; %white
fila3=[0,0,BlueMax]; %blue
J1=[fila1(1);fila2(1);fila3(1)];
J2=[fila1(2);fila2(2);fila3(2)];
J3=[fila1(3);fila2(3);fila3(3)];
J1I=interp1(x,J1,xi)';
J2I=interp1(x,J2,xi)';
J3I=interp1(x,J3,xi)';
jetflag=[J1I,J2I,J3I]; %(48,3)
%........................
if RedMax==1 & BlueMax==1
    %ANADO LOS EXTREMOS ROJO-FUERTE Y AZUL-FUERTE (los cojo del "jet" colormap):
    pal=[0.9375,0.8750,0.8125,0.7500,0.6875,0.6250,0.5625,0.5000]';
    palinv=flipud(pal);
    ceros=zeros(8,1);
    MAProjo=[palinv,ceros,ceros]; %(8 x 3)
    MAPazul=[ceros,ceros,pal]; %(8 x 3)
    jetflag=[MAProjo;jetflag;MAPazul]; %(64 x 3) donde "negativo-rojo, positivo-azul".
end
%........................
%OBTENGO TAMBIEN LA PALETA DE COLORES INVERSA:
jetflaginv=flipud(jetflag); %(64 x 3) donde "negativo-azul, positivo-rojo".
%........................
%............
% $$$ X=[-10:10];
% $$$ figure(1),plot(X),colormap(jetflag),colorbar('horiz')
% $$$ figure(2),plot(X),colormap(jetflaginv),colorbar('horiz')
%............
