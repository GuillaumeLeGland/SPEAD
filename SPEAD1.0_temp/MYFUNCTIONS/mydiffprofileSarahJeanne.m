%----------------------------------------------------------------------------
% $$$ function [Xout,varargout]=mydiffprofile(ndays,zmax,dz,mld,varname,ysurface)
%----------------------------------------------------------------------------

%............................................................................
ndays = 360; %[days]
zmax = 200; %[m]
dz = 2; %[m]
%............................................................................
T = ndays;
tau = T;
tmin = 1;
tmax = ndays;
dt = 1;
t = [dt:dt:tmax];
tt = t-tau;
wrad = (2*pi/T);
%............................................................................
mldmin = 10; %[m]
mldmax = 150; %[m]
A0 = mldmin;
A = (mldmax-mldmin)/2; %amplitud.
mld = A0 + A*( sin(wrad*tt) + 1); %MLD (seasonal)
%............................................................................
skzmin = 100; %[m2*day-1]
skzmax = 300; %[m2*day-1]
A0 = skzmin;
A = (skzmax-skzmin)/2; %amplitud.
skz = A0 + A*( sin(wrad*tt) + 1); %Surface kz (seasonal)
%............................................................................
varname = 'KZ';
%............................................................................
% $$$ mld = mldmax*ones(1,ndays); %For constant values (no seasonality)
% $$$ skz = skzmax*ones(1,ndays);
%............................................................................
ysurface = skz;
%............................................................................
%----------------------------------------------------------------------------

%*********************************
%MYDIFFPROFILE.m: Esta programa obtiene una matriz 2D (depth,time) 
%con perfiles de KZ (diff. turbulenta) o ST (sea temperature). Use:
%
%a) [KZ]=mydiffprofile(ndays,zmax,dz,mld,'KZ',skz); %for KZ
%
%b) [ST]=mydiffprofile(ndays,zmax,dz,mld,'ST',sst); %for ST.
%
%where:
%
% --------------------------------
% ndays: Temporal scale (ie. 365 days)
% zmax: Vertical scale (ie. 200 m)
% dz: Vertical resolution (ie. 5 m)
% mld: Vector with the mixed layer depth values (one per day).
% varname: 'KZ' (turbulence) or 'ST' (sea temperature)
% ysurface: Surface values (365 days) of "kz" or "st".
% --------------------------------
%*********************************
z=[0:dz:zmax-dz]';

%%%%%%%%%%%%%%%%%%%%
%DEFINE MIN AND MAX:
%%%%%%%%%%%%%%%%%%%%
if strcmp(varname,'ST')
    ymin=min(ysurface)-0.25; %scalar - global sea temperature minimum [C].
elseif strcmp(varname,'KZ')
    ymin=1; %scalar Global turbulent diffusion minimum [m2*day-1].
end
xmax=ysurface; %row-vect(1,ndays) daily maximum diffusivity [m2*day-1] or temp [C].
xmin=ymin*ones(1,ndays); %row-vect(1,ndays) daily minimum diffusivity [m2*day-1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE STEP OF THE PYCNOCLINE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------
%NOTA: Si defino xstepstar = xstep * (xmax - xmin), donde xmin y
%xmax son valores normalizados (entre 0 y 1), la relacion entre xstepstar
%y deltaZ de la pycnocline es:
%
%deltaZ = zmax * xstepstar^(-1);
%
%(see "egDiffProfile4.m")
%------------------------------
%....................
xstep=+20; %Original.
%....................
xstep=xstep*ones(1,ndays); %row-vect.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZE THE VARIABLES (VALUES BETWEEN 0 AND 1):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZE MLD AND DEPTHS BY THE MAXIMUM DEPTH OF THE DOMAIN:
mmld=mld/zmax; %row-vect.
zz=z/zmax; %col-vect.
zzmax=zmax/zmax; %scalar (equal to 1).
%.......................
ddz=dz/zmax; %scalar.
zzbis=[ddz:ddz:zzmax]';
[zz,zzbis]; %deben ser iguales.
%.......................
%NORMALIZE KZ BY ITS MAXIMUM VALUE AT THE SURFACE:
xxmax=xmax./xmax; %row-vect of maximum(s) KZ [adim]
xxmin=xmin./xmax; %row-vect of minimum(s) KZ [adim]
%.........

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET VERTICAL PROFILES FOR EACH DAY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xadims=[];
Xunits=[];
ExtraDepth=[];
Step=[];
for i=1:ndays
    %..............
    xmaxi=xmax(i); %absolute units [m2*day-1]
    xxmaxi=xxmax(i); %normalized values [0 - 1]
    xxmini=xxmin(i); %normalized values [0 - 1]
    xstepi=xstep(i); 
    %..............

    %====================================================
    %OBTENGO EL DELTA ENTRE VALOR MIN Y MAX NORMALIZADOS:
    %====================================================
    deltaxxi=xxmaxi-xxmini; %(between 0 and 1)

    %======================================
    %STANDARIZO THE STEP OF THE PYCNOCLINE:
    %======================================
    %---------------------
    %NOTE:
    %slope of pynocline = deltaZ / deltaX:
    %a) if deltaZ of pynocline is constant => slope decreases (more horizontal) when deltaX increase.
    %b) if deltaZ of pynocline increases with deltaX of pycnocline => slope increases (more vertical) when deltaX increase.
    %---------------------
    %................
% $$$     xxstepi=xstepi; %deltaZ of pynocline constant: slope decreases (more horizontal) when deltaX increase.
    %................
    xxstepi=(xstepi*deltaxxi); %deltaZ of pynocline increases with deltaX of pycnocline: slope increases (more vertical) when deltaX increase.
    %................
    if deltaxxi==0
	xxstepi=eps;
    end
    Step=[Step,xxstepi];
    
    %=================================================
    %GET DEPTH AT WICH I WANT THE PYCNOCLINE TO START:
    %=================================================
    %----------------------------------------------------------------------------
    %deltaZ = zmax * (xstepi * deltaxx)^(-1) = zmax * xxstep^(-1) = zmax/xxstep.
    %----------------------------------------------------------------------------
    %.............
    deltaZi = zmax./xxstepi; %delta Z(m) of pycnocline.
    extradepth=deltaZi/2;
    eextradepth=extradepth/zmax;
    %.............
    ExtraDepth=[ExtraDepth,extradepth];

    %====
    %MLD:
    %====
    if strcmp(varname,'ST')
	%..............
% $$$ 	mmldi=mmld(i);
	%..............
	mmldi=mmld(i)+eextradepth; 
	%..............
	hmld(i)=mmldi*zmax; %Remove normalization.
	%..............
	mldsst(i)=hmld(i); %save.
	%..............
    elseif strcmp(varname,'KZ')
	%..............
% $$$ 	mldssti=mldsst(i); %mld from SST + extradepht.
% $$$ 	mmldssti=mldssti/zmax;
% $$$ 	mmldi=mmldssti
	%..............
% $$$ 	mmldi=mmld(i);
	%..............
	mmldi=mmld(i)+eextradepth; 
	%..............
	hmld(i)=mmldi*zmax; %Remove normalization.
	%..............
	mldskz(i)=hmld(i); %save.
	%..............
    end
   
    %===================
    %SIGMOIDAL EQUATION:
    %===================
    %---------------------
    %NOTE:
    %slope of pynocline = deltaZ / deltaX:
    %a) beta1: if deltaZ of pynocline increases with deltaX of pycnocline => slope decreases (more horizontal) when deltaX increase.
    %b) beta2: => if deltaZ of pynocline is constant => slope decreases (more horizontal) when deltaX increase.
    %---------------------
    alfa = 2*(zz-mmldi)/zzmax;
    beta = exp(xxstepi*alfa);
    numer = xxmaxi+(xxmini*beta);
    denom = 1+beta;
    xx = numer./denom;
    x = xx*xmaxi; %recupero las dimensiones originales (remove normalization).
    
    %=========
    %STOCKAGE:
    %=========
    Xadims=[Xadims,xx(:)];
    Xunits=[Xunits,x(:)];

    %===============
    %VERTICAL PLOTS:
    %===============
% $$$     if i>200
% $$$     figure(1)
% $$$     plot(x,z*(-1),'-r.')
% $$$     axis([min(xsst) max(xsst), zmax*(-1) 0])
% $$$     grid on
% $$$     hold on
% $$$     pause(1)
% $$$     end
end
%ExtraDepth,pause

%???????????????????
% $$$ zmax
% $$$ deltaZi
% $$$ extradepth
% $$$ eextradepth
% $$$ subplot(2,2,1)
% $$$ plot(mmld)
% $$$ subplot(2,2,2)
% $$$ plot(h)
% $$$ pause
%???????????????????
    
%%%%%%%%
%PLOTS:
%%%%%%%%
%.....................................
XXadims=Xadims;
XXadims(1,1)=0;
XXadims(1,2)=1;
figure(2)
%.....................................
subplot(2,2,1)
imagesc(XXadims) %Normalized values [ 0 - 1 ]
colorbar('vertic')
% $$$ hold on
% $$$ hp=plot([1:ndays],mld/dz,'-',[1:ndays],hmld/dz,':');
% $$$ set(hp,'Color',[1 1 1],'LineWidth',[2]);
% $$$ hold off
xlabel('time [days]')
ylabel('depth [m]')
title('normalized units [ 0 - 1 ] n.d.')
%.....................................
subplot(2,2,2)
imagesc(Xunits) %Absolute values [0 - kzmax]
colorbar('vertic')
hold on
hp=plot([1:ndays],mld/dz,'-',[1:ndays],hmld/dz,':');
set(hp,'Color',[1 1 1],'LineWidth',[2]);
hold off
xlabel('time [days]')
ylabel('depth [m]')
title('absolute units [ 0 - kzmax ] m2*day-1')
pause(1)
%.....................................

%%%%%%%%
%OUTPUT:
%%%%%%%%
Xout=X;
varargout{1}=hmld;



