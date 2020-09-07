function [TDAYdaily,TDAYmonthly]=corrGLOBAL_daylength()

format short g
meses=['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

dt=365/12;
dlat=10;
xtick=[0:dt:365];xtick(1)=1;
ytick=[0:dlat:180];ytick(1)=1;
LATRG=[+90:-dlat:-90];

%*********************************************
%Programm: MYDAYLENGTH.m
%Este programa calcula la duracion del fotoperiodo en funcion de la
%latidud y dia del ano.
%Las equaciones las he obtenido de la pagina web:
%<http://www.gardenwithinsight.com/help100/00000254.htm> %para "Solar DayLength".
%<http://en.wikipedia.org/wiki/Declination#Sun> %para la "Solar Declination".
%
%Use: [TDAYdaily,TDAYmonthly]=mydaylength
%
%Outputs:
%
%TDAYdaily(180,365): daylength daily by latitude (fila1 es 90N, fila180 es 90S)
%TDAYmonthly(180,360,12): daylength monthly-averaged clim.
%*********************************************

%FACTOR DE CONVERSION DE LATITUD EN GRADOS A LATITUD EN RADIANTES:
fc=(2*pi)/360; %factor conversion grados2radianes.
LAT=[[-90:-1],[+1:+90]];
%LAT=[-89:+90];
BETA=fc*LAT; %latitud (en radianes).

%SOLAR DECLINATION:
DAY=[1:365];
TETA=-(fc*23.45)*cos(((2*pi)/365)*(DAY+10)); %solar declination (NOTA: Lleva "+10" para que el 21 Jun sea el dia mas largo del ano).

%DAYLENGTH (TDAY) PARA CADA LATITUD Y DIA DEL ANO:
TDAY=[];
for lati=BETA
    latitud=lati/fc;
    TDAYi=(24/pi)*acos(-tan(TETA)*tan(lati)); %daylenght para dayi.
    TDAY=[TDAY;TDAYi]; %(180,365)=(lat,day);
end

%æææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
%COJO SOLO LA PARTE "REAL" (A VECES SALEN NUMEROS COMPLEJOS, NO SE PQ!):
%æææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
TDAY=real(TDAY);

%%%%%%%
%PLOTS:
%%%%%%%
%NOTA: Comparar "figure(1)" con "/home/svallina/SER/FOTOS/daylength.GIF"
%(que viene de la web <http://www.physicalgeography.net/fundamentals/6i.html>)
figure(1)
subplot(1,2,1)
plot(DAY,TDAY(160,:),'.b-') %70N
hold on
plot(DAY,TDAY(150,:),'.g-') %60N
hold on
plot(DAY,TDAY(140,:),'.r-') %50N
hold on
plot(DAY,TDAY(120,:),'.y-') %30N
hold on
plot(DAY,TDAY(91,:),'.m-') %1N
hold off
grid on
axis([0 365, 0 24])

subplot(1,2,2)
plot(DAY,TDAY(160,:),'.b-') %70N
hold on
plot(DAY,TDAY(150,:),'.g-') %60N
hold on
plot(DAY,TDAY(140,:),'.r-') %50N
hold on
plot(DAY,TDAY(120,:),'.y-') %30N
hold on
plot(DAY,TDAY(91,:),'.m-') %1N
hold off
grid on
set(gca,'Ylim',[0 24],'YTick',[0:4:24],'YTickLabel',[0:4:24],'Xlim',[1 365],'XTick',xtick,'XTickLabel',meses,'FontSize',[4])
legend('70N','60N','50N','30N','1N')
%OUTPUT:
TDAYdaily=flipud(TDAY); %fila1=90N; fila180=90S.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONSTRUYO MATRIZ 3D DE MONTHLY MEANS (180,360,12):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mlims=[0,31,28,31,30,31,30,31,31,30,31,30,31];
for k=1:12
    mlim0=sum(Mlims(1:k))+1;
    mlim1=sum(Mlims(1:k+1));
    [mlim0,mlim1];
    TDAYk=snanmean(TDAY(:,mlim0:mlim1),2); %empieza por el HS (la primera pos es 90S).
    TDAYk=flipud(TDAYk); %empiezo por el HN (primera pos = 90N).
    TDAYmonthly(:,:,k)=TDAYk*ones(1,360);
end

figure(2)
for k=1:12
    subplot(3,4,k)
    TDAYmonthlyk=TDAYmonthly(:,:,k);
    TDAYmonthlyk(1,1)=0;
    TDAYmonthlyk(1,2)=24;
    imagesc(TDAYmonthlyk)
    hc=colorbar('horiz');
    set(gca,'Ytick',ytick,'Yticklabel',LATRG,'Fontsize',[6])
    set(hc,'Xlim',[0 24],'Xtick',[0:4:24],'Xticklabel',[0:4:24],'Fontsize',[6])
end
pause(1)
close(1)
close(2)

%%%%%%
%SAVE:
%%%%%%
% $$$ TDAY=TDAYmonthly; %cambio el nombre pq es mas corto (OJO! no confundir la
% $$$                   %nueva variable TDAY (que ahora es 3D monthly) con la
% $$$                   %original TDAY(180,365)(lat,days).
% $$$ save /SERVAL/SER24/PROGRAMMING/MATLAB/RESULTADOS/CORRELACION/GLOBAL/VARIABLES/DAYLENGTH/TDAY.mat TDAY
