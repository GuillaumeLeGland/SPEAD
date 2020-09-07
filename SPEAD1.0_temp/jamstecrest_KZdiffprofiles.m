%function [Xout] = jamstecrest_KZdiffprofiles(ndays,zmax,dz,mld,ymax)
% Now takes a deep value into account (Le Gland, 24/10/2019)
function [Xout] = jamstecrest_KZdiffprofiles(ndays,zmax,dz,mld,ymax,ydeep)

%z = [0:dz:zmax-dz]; %row vector.
z = [0:dz:zmax]; % Depths at each box boundaries (Le Gland, 06/09/2019)
z = z(:); %make it column vector.

%%%%%%%%%%%%%%%%%%%%
%DEFINE MIN AND MAX:
%%%%%%%%%%%%%%%%%%%%
xmax = ymax; %row-vect(1,ndays) daily maximum diffusivity [m2*day-1]
%xmin = ones(1,ndays); %row-vect(1,ndays) daily minimum diffusivity [m2*day-1]
xmin = ydeep*ones(1,ndays);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE STEP OF THE PYCNOCLINE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xstep=+20;
xstep=xstep*ones(1,ndays); %row-vect.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZE THE VARIABLES (VALUES BETWEEN 0 AND 1):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmld=mld/zmax; %row-vect.
zzmax=zmax/zmax; %scalar.
zz=z/zmax; %col-vect.
xxmax=xmax./xmax; %row-vect. of maximums SST [adim]
xxmin=xmin./xmax; %row-vect. of minimums SST [adim]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET VERTICAL PROFILES FOR EACH DAY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[];
Xadim=[];
ExtraDepth=[];
Step=[];
for i = 1:ndays
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
	      xxstepi=eps; %almost zero (to avoid infinite values when dividing by xxstepi).
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
    mmldi=mmld(i)+eextradepth; 
    h(i)=mmldi*zmax;
    
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
end

%%%%%%%%
%OUTPUT:
%%%%%%%%
Xout = X;



return

