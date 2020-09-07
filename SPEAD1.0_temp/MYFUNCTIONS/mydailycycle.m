function [PDF]=mydailycycle(sunrise,dt)

%*************************************************************************
%Programa: MYDAILYCYCLE.m
%Este programa calcula el ciclo solar de radiacion entre sunrise y sunset
%usando una funcion sinusoidal. Da valores normalizados (entre 0 y 1).
%
%Use: [PDF]=mydailycycle(sunrise,dt)
%
% sunrise: hora a la que aparce el sol [h].
% dt: time step para calcular la curva del ciclo solar [h].
%*************************************************************************

daytime=[0:dt:24]; %[h]
n=length(daytime);

delta=6-sunrise;
sunset=18+delta;
%[sunrise,sunset]

daylength=(sunset-sunrise);
t0=sunset-daylength;

PDF=sin(pi/daylength*(daytime-t0)); 

I=find(daytime<sunrise);
J=find(daytime>sunset);
PDF(I)=0;
PDF(J)=0;
%..................................
% $$$ figure(10)
% $$$ plot(daytime,PDF,daytime,PDF,'r.')
% $$$ axis([0 24, 0 1])
% $$$ set(gca,'Xtick',[0:1:24])
% $$$ grid on
%..................................
