function []=myspatialcorrmaps(VARX,VARY,h,MAP,fignum)
load /home/svallina/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/CORRELACION2/GLOBAL/Land.txt
%NOTA: Si uso "mensual" solo funciona bien para h=3 (window 3 x 3), no se pq.

[Ex,Ey,Exy,Vx,Vy,N]=myspatialcorrmaps_estadisticos(VARX,VARY,h,Land,'mensual');
%[Ex,Ey,Exy,Vx,Vy,N]=myspatialcorrmaps_estadisticos(VARX,VARY,h,Land,'anual');
[RP,RPsig,TP,LtotP]=myspatialcorrmaps_rpearson(Ex,Ey,Exy,Vx,Vy,N,Land,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT DEL MAPA DE CORRELACIONES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myglobalmap(RP,'R',Land,jet,fignum);
