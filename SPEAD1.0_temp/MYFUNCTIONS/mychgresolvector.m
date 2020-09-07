function [xave] = mychgresolvector(xdata,wj)

%*********************************************************************
%MYCHGRESOL.m: Este programa reduce la resolucion de x(m) a xout(n), 
%donde "m > n" (o sea, disminuyo la resolucion: paso de muchos pixeles
%pequenos a pocos pixeles grandes). Se basa en promediar sobre grupos de
%pixeles pequenos para obtener uno grande para cada grupo, de modo que la
%size de x y xout debe ser un multiplo entero.
%
%"wj" es el numero de pixels que entran la window.
%
% [xave] = mychgresolvector(xdata,wj)
%
%*********************************************************************
msize = length(xdata);
counter = 0;
for j = 1:wj:msize 
    counter = counter + 1; 
    window = [j:j+(wj-1)];
    xwin = xdata(window);
    xave(counter) = nanmean(xwin);
end
return

