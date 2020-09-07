function [Xout,varargout]=mydiffprofile(ndays,zmax,dz,mld,varname,ysurface)
global mldsst 

%*********************************
%MYDIFFPROFILE.m: Esta programa obtiene una matriz 2D (depth,time) con
%perfiles de KZ (diff. turbulenta) o ST (sea temperature).
%Use:
%
%a) mydiffprofile(zmax,dz,mld,'KZ'); %for KZ
%
%b) mydiffprofile(zmax,dz,mld,'ST',sst); %for ST.
%
%ysurface: Surface values (365 days)
%*********************************
%.................
z=[0:dz:zmax-dz]'; %USAR ESTE!!
%.................
% $$$ z=[dz:dz:zmax]';
%.................

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
mmld=mld/zmax; %row-vect.
zzmax=zmax/zmax; %scalar.
zz=z/zmax; %col-vect.
%.......................
ddz=dz/zmax; %scalar.
zzbis=[ddz:ddz:zzmax]';
[zz,zzbis]; %deben ser iguales.
%.......................
xxmax=xmax./xmax; %row-vect of maximums SST [adim]
xxmin=xmin./xmax; %row-vect of minimums SST [adim]
%.........

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET VERTICAL PROFILES FOR EACH DAY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[];
Xadim=[];
ExtraDepth=[];
Step=[];
for i=1:ndays
    %..............
    xmaxi=xmax(i);
    xxmaxi=xxmax(i);
    xxmini=xxmin(i);
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
	h(i)=mmldi*zmax;
	mldsst(i)=mmldi*zmax; %save.
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
	h(i)=mmldi*zmax;
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
    alfa=2*(zz-mmldi)/zzmax;
    beta=exp(xxstepi*alfa);
    numer=xxmaxi+xxmini*beta;
    denom=1+beta;
    xx=numer./denom;
    x=xx*xmaxi; %recupero las dimensiones originales.
    
    %=========
    %STOCKAGE:
    %=========
    Xadim=[Xadim,xx(:)];
    X=[X,x(:)];

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
% $$$ XXadim=Xadim;
% $$$ XXadim(1,1)=0;
% $$$ XXadim(1,2)=1;
% $$$ 
% $$$ figure(2)
% $$$ subplot(2,2,1)
% $$$ imagesc(XXadim)
% $$$ colorbar('horiz')
% $$$ hold on
% $$$ hp=plot([1:ndays],mld/dz,'-',[1:ndays],h/dz,':');
% $$$ set(hp,'Color',[1 1 1],'LineWidth',[2]);
% $$$ hold off
% $$$ 
% $$$ subplot(2,2,2)
% $$$ imagesc(X)
% $$$ colorbar('horiz')
% $$$ hold on
% $$$ hp=plot([1:ndays],mld/dz,'-',[1:ndays],h/dz,':');
% $$$ set(hp,'Color',[1 1 1],'LineWidth',[2]);
% $$$ hold off
% $$$ pause
%.....................................

%%%%%%%%
%OUTPUT:
%%%%%%%%
Xout=X;
varargout{1}=h;



