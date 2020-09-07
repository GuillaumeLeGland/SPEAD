function [RS,RSsig]=mycorrmaps(VARX,VARY,varargin)
%function []=mycorrmaps(VARX,VARY,MAP,fignum)
%**********************************************************
%PROGRAMA: CORR.m (17 May 2004. USAR ESTE!)
%
%Este programa calcula el mapa de correlacion 2-D entre
%dos matrices 3D (180,360,12), obtenido pixel-a-pixel 
%para series temporales de 12 meses (donde los puntos son el valor
%promedio de una ventana deslizante 7x7).
% 
%Uso: 
% a) [RS,RSsig]=mycorrmaps(VARX,VARY,MAP,fignum) %grafico mapa de corrs.
% b) [RS,RSsig]=mycorrmaps(VARX,VARY) %NO grafico mapa de corrs.
%
%donde:
%VARX(180,360,12).
%VARY(180,360,12).
%MAP = colormap (64,3). (ie. "jet").
%fignum = figure number.
%**********************************************************
%.............................................
% $$$ load /home/svallina/SERVAL/SER24/DATA/GLOBAL/DMS/DMSc/KETTLE/DMS_KETTLE_2000.mat
% $$$ load /home/svallina/SERVAL/SER24/DATA/GLOBAL/CCN/MODIS/CCNmodclim0204.mat
% $$$ VARX=DMS_KETTLE_2000;
% $$$ VARY=CCNmodclim0204;
%.............................................

% $$$ addpath /home/svallina/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/STOOLBOX
% $$$ addpath /home/svallina/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYFUNCTIONS/
% $$$ path=genpath('/home/svallina/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/TOOLBOX/');
% $$$ addpath(path)
% $$$ load /home/svallina/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/CORRELACION2/GLOBAL/Land.txt
format short g

[GLOBE,Land]=myglobeland;

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EL NUMERO DE VARARGIN DEFINE SI GRAFICO EL MAPA DE CORRELACIONES O SOLO
%ME QUEDO CON LA MATRIZ 2D(lat,long) DE CORRELACIONES (Y LA GRAFICO LUEGO YO):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numvar=length(varargin);
if numvar==2 %grafico mapa de correlaciones.
    MAP=varargin{1};
    fignum=varargin{2};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRETRATAMIENTO DE LAS MATRICES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n,p]=size(VARX);
%============================
%APLICO LA MASCARA DE "LAND":
%============================
for k=1:p
    Xk=VARX(:,:,k);
    Yk=VARY(:,:,k);
    Xk(Land)=nan;
    Yk(Land)=nan;
    VARX(:,:,k)=Xk;
    VARY(:,:,k)=Yk;
end

%===============================================
%BUSCO LAS POSICIONES CON NO-DATA Y LES DOY NAN:
%===============================================
Ineg=find(VARX<0);
Jneg=find(VARY<0);
VARX(Ineg)=nan;
VARY(Jneg)=nan;
clear Ineg Jneg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GRAFICO LAS CARTAS MENSUALES DE VARX e VARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ mymonthlymaps(VARX,'VARX',Land,jet,1);
% $$$ mymonthlymaps(VARY,'VARY',Land,jet,2);
% $$$ pause(1)
% $$$ close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PONGO NaN EN LAS POS(i,j) DONDE VARX o VARY TIENEN NAN:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F=find(isnan(VARX)==1 | isnan(VARY)==1);
VARX(F)=nan;
VARY(F)=nan;
clear F

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO LOS MAPAS DE CORRELACION 2D (SPEARMAN):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================================
%COEF. CORRELACION TEMPORAL (SPEARMAN):
%======================================
%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
%OBTENGO PROMEDIOS MENSUALES PARA CADA PIXEL (PROMEDIO EN EL ESPACIO
%EN BASE A LA WINDOW PARA CADA MES):
%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
w=3; %w=3 para window 7x7
%[VARXm,VARXstd]=srunmean(VARX,w);
%[VARYm,VARYstd]=srunmean(VARY,w);
[VARXm,VARXstd]=myrunmean(VARX,w,'Espacial');
[VARYm,VARYstd]=myrunmean(VARY,w,'Espacial');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PONGO NaN EN LAS POS(i,j) DONDE VARXm o VARYm TIENEN NAN:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F=find(isnan(VARXm)==1 | isnan(VARYm)==1);
VARXm(F)=nan;
VARYm(F)=nan;
clear F

%æææææææææææææææææææææææææææææææææææææææææææ
%CALCULO LA CORRELACION TEMPORAL (SPEARMAN):
%æææææææææææææææææææææææææææææææææææææææææææ
[RS,RSsig,TS]=mycorrmaps_rspearman(VARXm,VARYm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT DEL MAPA DE CORRELACIONES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numvar==2
    myglobalmap(RS,'RS',Land,MAP,fignum);
    %myglobalmap(RSsig,'RS',Land,MAP,fignum);
end

%********************
toc
