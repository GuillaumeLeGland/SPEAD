function [jetflag]=jetflagcolorbar()
close all
clear all

x=[1,32,64];
xi=[1:64];
fila1=[1,0,0]; %red
fila2=[1,1,1]; %white
fila3=[0,0,1]; %blue
J1=[fila1(1);fila2(1);fila3(1)];
J2=[fila1(2);fila2(2);fila3(2)];
J3=[fila1(3);fila2(3);fila3(3)];
J1I=interp1(x,J1,xi)';
J2I=interp1(x,J2,xi)';
J3I=interp1(x,J3,xi)';

jetflag=[J1I,J2I,J3I]

%............
% $$$ X=[-10:10];
% $$$ plot(X),colormap(jetflag),colorbar('horiz')
%............
