function [EFLUX,EDAY]=myeday()

tic
%***************************************************************************
%Calcula la 'Solar Irrandiance Flux (E)' para cada latitud en funcion del
%ciclo de rotacion de la Tierra (el primer dia del ciclo es el equinocio
%de boreal-autom (21 sept)):
%***************************************************************************

%MODEL KEYS:
keyTimeScaleAve = 'Weekly';
% $$$ keyTimeScaleAve = 'Monthly';

if strcmp(keyTimeScaleAve,'Monthly')
    psize = 12; %Ave 30 days.
    Mlims = [0,31,28,31,30,31,30,31,31,30,31,30,31];
elseif strcmp(keyTimeScaleAve,'Weekly')
    psize = 72; %Ave 05 days.
    MlimsJan = [5,5,5,5,5,6];
    MlimsFeb = [5,5,5,5,5,3];
    MlimsMar = [5,5,5,5,5,6];
    MlimsApr = [5,5,5,5,5,5];
    MlimsMay = [5,5,5,5,5,6];
    MlimsJun = [5,5,5,5,5,5];
    MlimsJul = [5,5,5,5,5,6];
    MlimsAug = [5,5,5,5,5,6];
    MlimsSep = [5,5,5,5,5,5];
    MlimsOct = [5,5,5,5,5,6];
    MlimsNov = [5,5,5,5,5,5];
    MlimsDec = [5,5,5,5,5,6];
    Mlims = [0,MlimsJan,MlimsFeb,MlimsMar,MlimsApr,MlimsMay,MlimsJun,MlimsJul,MlimsAug,MlimsSep,MlimsOct,MlimsNov,MlimsDec];
end

%FACTOR DE CONVERSION DE LATITUD EN GRADOS A LATITUD EN RADIANTES:
lonmax = 180;
latmax = 360;
fc=(2*pi)/latmax; %factor conversion grados2radianes.

%DEFINO EL CICLO DE TRANSLACION DE LA TIERRA (1-365) COMO RADIANES (0-2*pi):
drot=(2*pi)/365;
ROT=[0:drot:(2*pi)-drot];

%CALCULO QUE LATITUD (FI*) ESTA RECIBIENDO LOS RAYOS DEL SOL DE FORMA
%PERPENDICULAR (VARIA SEGUN EL MOMENTO DEL ANO):
alfa=23.45; %angulo entre el plano del Sol y el plano de translacion de la Tierra ("Declination angle").
ALFAstar=sin(ROT)*alfa; %angulo entre el plano del Sol y el equador de la Tierra.
BETA=90+ALFAstar; %angulo entre el plano del sol y el eje de rotacion de la Tierra.
FIstar=90-BETA;

%CALCULO EL SOLAR ELEVATION ANGLE PARA CADA LATITUD A CADA MOMENTO DEL ANO:
%------------------------------------------------------------
%NOTA: BETA(1)=0, corresponde a 21 Sept (equinocio de OTONO).
%------------------------------------------------------------
FI=[+90:-1:-89]'; %latitudes.
%FI=[+90:-0.25:-90]';
m=length(FI);
PSI=[];
for j=1:365 %i=1 corresponde a 21 Sept (equinocio de OTONO).
    %..................................................
% $$$     betaj=BETA(j);
% $$$     PSIj=betaj-FI; %ELEVATION ANGLES PARA EL DIA j.
    %..................................................
    alfastarj=ALFAstar(j);
    PSIj=90-(alfastarj+FI); %ELEVATION ANGLES PARA EL DIA j.
    %..................................................
    PSI=[PSI,PSIj];
end

%CALCULO EL MAXIMUM SOLAR IRRADIANCE FLUX (E) (SE DEFINE COMO EL FLUJO DE
%RADIACION QUE RECIBE UNA SUPERFICIE 'PERPENDICULAR' A LOS RAYOS):
Ro=149.6; %[Gm = 10^9 m] distancia promedio entre el Sol y la Tierra.
So=1368; %[W*m-2] = solar constant.
R=[];
for j=1:365 %i=1 corresponde a 21 Sept (equinocio de OTONO).
    jj=(j+78); %days from 4 July (maximum distance between Sun and Earth).
    if jj>=365 %dia 4 july.
	jj=jj-365;
    end
    jj
    jjrad=jj*(2*pi/365); %days from 4 July (en radianes).
    Rj=Ro+2.5*cos(jjrad); %[m] distancia promedio entre el Sol y la Tierra para el dia j.
    R=[R,Rj];
end
E=So*(Ro./R).^2; %[W*m-2].
%.........
plot(R)
close
%.........

%CALCULO EL SOLAR IRRADIANCE FLUX (E*) A CADA LATITUD Y DIA DEL ANO:
Estar=[];
for j=1:365 %i=1 corresponde a 21 Sept (equinocio de OTONO).
    PSIj=PSI(:,j); %elevation angles para dia-j.
    Ej=E(j);
    Estarj=sin(PSIj*fc)*Ej;
    Estar=[Estar,Estarj];
end
I=find(Estar<0);
Estar(I)=nan; %zonas de sombra (donde no da el sol).
MAP=jet;
MAP(1,:)=[0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AVERAGED DAILY SOLAR FLUX:
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULO LA DAILY SOLAR IRRADIANCE:
%--------------------------------------
%NOTA:
%Eday = [(So/pi)*(Ro/R)^2]*[ho*sin(fi)*sin(del) + cos(fi)*cos(del)*sin(ho)] =
%     = [E/pi]*[ho*sin(fi)*sin(del) + cos(fi)*cos(del)*sin(ho)].
%
%where:
%
%E: Maximum solar irradiance.
%ho = arcos(-tan(fi)*tan(del)): hour angle in radians (fi y del deben estar en radianes).
%d : days from 21 June.
%drad = d*(2*pi/365): days from 21 June (in radians).
%del = [23.45*(2*pi/latmax)]*cos(drad): declination angle (in radians).
%.......................
%Example1:
% $$$ fi=49.25
% $$$ firad=fi*(2*pi/latmax);
% $$$ d=1 %21 June.
% $$$ drad=d*(2*pi/365)
% $$$ del=(23.45*fc)*cos(drad)
% $$$ ho=acos(-tan(firad)*tan(del))
% $$$ Ro=149.6;
% $$$ R=151.892;
% $$$ Edayj=((So/pi)*(Ro/R)^2)*(ho*sin(firad)*sin(del) + cos(firad)*cos(del)*sin(ho)) %= 485.74 [W*m-2]
%......................................
%Example2:
% $$$ fi=49.25
% $$$ firad=fi*(2*pi/latmax);
% $$$ d=1+182 %21 Dec.
% $$$ drad=d*(2*pi/365)
% $$$ del=(23.45*fc)*cos(drad)
% $$$ ho=acos(-tan(firad)*tan(del))
% $$$ Ro=149.6;
% $$$ R=151.892;
% $$$ Edayj=((So/pi)*(Ro/R)^2)*(ho*sin(firad)*sin(del) + cos(firad)*cos(del)*sin(ho)) %= ? [W*m-2]
%......................................
%Example3:
% $$$ fi=90
% $$$ firad=fi*(2*pi/latmax);
% $$$ d=1+182 %21 Dec.
% $$$ drad=d*(2*pi/365)
% $$$ del=(23.45*fc)*cos(drad)
% $$$ ho=acos(-tan(firad)*tan(del))
% $$$ Ro=149.6;
% $$$ R=151.892;
% $$$ Edayj=((So/pi)*(Ro/R)^2)*(ho*sin(firad)*sin(del) + cos(firad)*cos(del)*sin(ho)) %= ? [W*m-2]
%--------------------------------------
Eday=[];
D=[];
for j=1:365 %i=1 corresponde a 21 Sept (equinocio de OTONO).
    FIrad=FI*fc;
    ej=E(j); %Maximum solar irradiance para day-j.
    d=j+91; %from 21 Jun (1st day of NH summer).
    if d>=365
	d=d-365;
    end
    d
    drad=d*(2*pi/365); %days from 21 June (in radians).
    del=(23.45*fc)*cos(drad); %declination angle (in radians).
    Ho=acos(-tan(FIrad)*tan(del)); %hour angle (in radians).
    Edayj=(ej/pi)*(Ho.*sin(FIrad)*sin(del) + cos(FIrad).*cos(del).*sin(Ho));
    Eday=[Eday,Edayj(:)];
    D=[D;d];
end
Eday=real(Eday);
%...................................
%test: para D=0,
%a) at lat=50N (aprox lati=40) ==> Eday=485.74.
%b) at lat=Equator (aprox lati=90) ==> Eday<400.
[D,Eday(40,:)',Eday(90,:)',]; 
%...................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONSTRUYO MATRIZ 3D DE MONTHLY MEANS (180,360,12):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MAX. TOTAL FLUX:
Estar1Jan=[Estar(:,102:end),Estar(:,1:101)]; %hago empezar a Estar en el 1-Jan.
for k=1:psize
    mlim0=sum(Mlims(1:k))+1;
    mlim1=sum(Mlims(1:k+1));
    Efluxk=nanmean(Estar1Jan(:,mlim0:mlim1),2);
    EFLUX(:,:,k)=Efluxk*ones(1,latmax);
end

%AVERAGE DAILY INSOLATION:
Eday1Jan=[Eday(:,102:end),Eday(:,1:101)]; %hago empezar a Estar en el 1-Jan.
for k=1:psize
    mlim0=sum(Mlims(1:k))+1;
    mlim1=sum(Mlims(1:k+1));
    Edayk=nanmean(Eday1Jan(:,mlim0:mlim1),2);
    EDAY(:,:,k)=Edayk*ones(1,latmax);
end

%%%%%%%
%PLOTS:
%%%%%%%
meses=['J','F','M','A','M','J','J','A','S','O','N','D','J']';
meses2=['S','O','N','D','J','F','M','A','M','J','J','A','S']';
figure(1)
plot(FIstar) %latitud que recibe los rayos de forma perpendicular a cada momento del ano.

figure(2)
EEstar=Estar;
EEstar(1,1)=-30;
imagesc(EEstar)
LATRG=[+90:-10:-90];
dt=365/12;
dlat=10*floor(m/lonmax);
xtick=[0:dt:365];xtick(1)=1;
ytick=[0:dlat:m];ytick(1)=1;
set(gca,'Xtick',xtick,'Xticklabel',meses2,'Ytick',ytick,'Yticklabel',LATRG)
mycolorbar('horiz')
colormap(MAP)

figure(3)
EEday=Eday;
EEday(1,1)=-30;
imagesc(EEday)
LATRG=[+90:-10:-90];
dt=365/12;
dlat=10*floor(m/lonmax);
xtick=[0:dt:365];xtick(1)=1;
ytick=[0:dlat:m];ytick(1)=1;
set(gca,'Xtick',xtick,'Xticklabel',meses2,'Ytick',ytick,'Yticklabel',LATRG)
mycolorbar('horiz')
colormap(MAP)

%POR BANDAS DE LATITUD:
figure(4)
c=0;
for j=10:20:170
    c=c+1;
    subplot(3,3,c)
    %...................................................
% $$$     plot([1:365],Estar(j,:)) %empieza el 21 Sep
% $$$     axis([1 365, 0 1500])
% $$$     set(gca,'Xtick',xtick,'Xticklabel',meses2,'Fontsize',[6]);
    %...................................................
    plot([1:365],Estar1Jan(j,:)) %empieza el 1 Jan
    axis([1 365, 0 1500])
    set(gca,'Xtick',xtick,'Xticklabel',meses,'Fontsize',[6]);
    %...................................................
    lat=(j-90)*(-1);
    title(['LATITUDE: ',num2str(lat)])
end

figure(5)
c=0;
for j=10:20:170
    c=c+1;
    subplot(3,3,c)
    %...................................................
% $$$     plot([1:365],Eday(j,:)) %empieza el 21 Sep
% $$$     axis([1 365, 0 1500])
% $$$     set(gca,'Xtick',xtick,'Xticklabel',meses2,'Fontsize',[6]);
    %...................................................
    plot([1:365],Eday1Jan(j,:)) %empieza el 1 Jan
    axis([1 365, 0 1500])
    set(gca,'Xtick',xtick,'Xticklabel',meses,'Fontsize',[6]);
    %...................................................
    lat=(j-90)*(-1);
    title(['LATITUDE: ',num2str(lat)])
end

if psize == 12
    nrows = 3;
    ncols = 4;
elseif psize == 72
    nrows = 3 + 5;
    ncols = 4 + 5;
end

figure(6)
for k=1:psize
    subplot(nrows,ncols,k)
    EFLUXk=EFLUX(:,:,k);
    EFLUXk(1,1)=0;
    EFLUXk(1,2)=1500;
    imagesc(EFLUXk)
    mycolorbar('horiz')
    set(gca,'Ytick',ytick,'Yticklabel',LATRG,'Fontsize',[6])
end

figure(7)
for k=1:psize
    subplot(nrows,ncols,k)
    EDAYk=EDAY(:,:,k);
    EDAYk(1,1)=0;
    EDAYk(1,2)=560;
    imagesc(EDAYk)
    mycolorbar('horiz')
    set(gca,'Ytick',ytick,'Yticklabel',LATRG,'Fontsize',[6])
end

%%%%%%%%%%%
%MASK LAND:
%%%%%%%%%%%
[GLOBE,Land]=myglobeland;
for k=1:psize
    %......................
    EDAYk=EDAY(:,:,k);
    %......................
    EDAYk(Land)=0;
    %......................
    EDAYocean(:,:,k)=EDAYk;
    %......................
end

%%%%%%
%SAVE:
%%%%%%
%.........................................................
% $$$ save /home/svallina//SERVAL/SER24/DATA/GLOBAL/LIGHT/EDAY/EDAYtoa.mat EDAYtoa
%.........................................................
% $$$ save Eflux.txt -ascii Estar %empieza el 21 Sept.
% $$$ save Eday.txt -ascii Eday %empieza el 21 Sept.
% $$$ save EFLUX.mat EFLUX % (matriz 3D) empieza el 1 Jan.
% $$$ save EDAY.mat EDAY % (matriz 3D) empieza el 1 Jan.
%.........................................................
%***********************************
toc
display('end: myeday.m')
disp('pause')
pause
close all
