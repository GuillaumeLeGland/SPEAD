function [W]=mySinkingProfile(ndays,zmax,dz,mld,wmin,wmax,varargin)

%*********************************
%MYSINKINGPROFILE.m: Esta programa obtiene una matriz 2D (depth,time) con
%perfiles de W (vertical sinking) en funcion del gradiente de densidad
%(en la picnoclina se reduce el sinking).
%
%Use: [W]=mySinkingProfile(zmax,dz,mld,wmin,wmax,varargin)
%
%*********************************
z=[0:dz:zmax-dz];

kzmax=ones(1,ndays);
kzmin=zeros(1,ndays);
kzgrad=-20; %picnocline gradient.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADIMENSIONALIZO LOS TERMINOS (FROM CROPP CODE, VALORES ADIM ENTRE O Y 1):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmld=mld/zmax;
zzmax=zmax/zmax;
zz=z/zmax;
%.......................
ddz=dz/zmax;
zzbis=[ddz:ddz:zzmax];
[zz,zzbis]; %deben ser iguales.
%.......................
kkzmax=kzmax./kzmax; %maximum diffusivity (adim)
kkzmin=kzmin./kzmax;   %minimum diffusivity (adim)
kkzgrad=kzgrad; %picnocline gradient.
KZadim=[];

for i=1:ndays
    mmldi=mmld(i);
    kzmaxi=kzmax(i);
    kkzmaxi=kkzmax(i);
    kkzmini=kkzmin(i);
    alfa=2*(zz-mmldi)/zzmax;
    beta=exp(kkzgrad*(kkzmini-kkzmaxi)*alfa);
    numer=kkzmaxi+kkzmini*beta;
    denom=1+beta;
    kkz=numer./denom;
    KZadim=[KZadim,kkz(:)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULO EL GRADIENTE VERTICAL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradKZ=diff(KZadim,1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZO (0-1) EL GRADIENTE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxgradKZ=max(abs(gradKZ(:)));
gradKZadim=gradKZ/maxgradKZ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO LAS VELOCIDADES DE SINKING:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmax=max(gradKZadim(:));
xmin=min(gradKZadim(:));

W = wmax - (wmax-wmin)*((xmax-gradKZadim)./(xmax-xmin));

%RECUPERO SIZE ORIGINAL:
W=[W;W(end,:)];

%%%%%%%%
%PLOTS:
%%%%%%%%
% $$$ figure(1)
% $$$ subplot(2,2,1)
% $$$ imagesc(KZadim)
% $$$ colorbar('horiz')
% $$$ 
% $$$ figure(1)
% $$$ subplot(2,2,2)
% $$$ imagesc(W)
% $$$ colorbar('horiz')
% $$$ 
% $$$ subplot(2,2,3)
% $$$ imagesc(gradKZ)
% $$$ colorbar('horiz')
% $$$ 
% $$$ subplot(2,2,4)
% $$$ imagesc(gradKZadim)
% $$$ colorbar('horiz')
% $$$ 
