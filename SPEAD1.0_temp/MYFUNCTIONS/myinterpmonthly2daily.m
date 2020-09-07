function [VARI]=bats_myinterpmonthly2daily(VAR);
VAR=VAR(:)'; %fila
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INTERPOLO LOS 12 VALORES DE VAR PARA OBTENER 365 VALORES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------------
%NOTA:
%Considero que los 12 valores de VAR estan centrados en el dia 15 de cada
%mes. De modo que tengo valores del 15 Jan hasta el 15 Dec (335 dias si
%interpolo). Voy a anadir el dato del 15 Dec del ano anterior y del 15
%Jan ano posterior para asi tener 395 dias interpolados. Cortando esta
%serie interpolada entre los dias 16 y 380 obtengo 365 dias interpolados
%(es decir del 1 Jan hasta el 31 Dec).
%------------------------------------------------------------------------
VVAR=[VAR(12),VAR,VAR(1)]; %anado Dec ano anterior y Jan ano siguiente.
dx=13/395; %14 meses (13 segmentos de 30 dias) dividido entre 335 del ano (del 15 Jan al 15 Dec) mas 60 dias anadidos (del 15 Dec ano anterior al 15 Jan ano posterior).
x=[1:14];
xI=[1:dx:14];
%VVARI=interp1(x,VVAR,xI,'spline');
VVARI=interp1(x,VVAR,xI,'cubic'); %para "bats_dataforSIMO".
VARI=VVARI(16:380); %quito del 15-31 Dec ano anterior y del 1-15 Jan ano
                    %posterior (asi me quedo entre 1 Jan y 31 Dec del ano actual).
figure(1)
plot([1:14],VVAR,'.b-')
axis([-inf +inf -inf +inf])
figure(2)
plot([1:396],VVARI,'.r-')
axis([-inf +inf -inf +inf])

figure(3)
plot([1:12],VAR,'.b-')
axis([-inf +inf -inf +inf])
figure(4)
plot([15:350],VARI(15:350),'.r-') %del 15 Jan al 15 Dec.
axis([-inf +inf -inf +inf])

%figure(5),plot([1:365],VARI*(-1),'.r-'),axis([0 365, -250 0])
figure(5),plot([1:365],VARI,'.r-'),axis([0 365, 0 +inf])
pause(1)

