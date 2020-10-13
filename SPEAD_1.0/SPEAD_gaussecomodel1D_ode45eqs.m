function [Vdot] = SPEAD_gaussecomodel1D_ode45eqs(iTime,V0)
global galfa gbeta 
global gzmax kgz mz betaz mpower 
global mp Isat InhFac numutx numuty
global alp0 mup0 
global amup aknp  
global Q10a Q10h
global ntot0 
global temp0 temp 
global jcounter  
global keyPhysics
global omePhy epsZoo omeZoo md 
global deltat
global zdepths ndepths deltaz
global jday 
global Ixave Iyave Ixxvar Iyyvar Ixycov
global Iphy Izoo Idin Ipon Ibox
global KZ
global parz0
global kw wsink
%................................................................................... 
global todedotday 
%................................................................................... 
global todedotout 
global UXYout GXYout 
%...................................................................................
global FPHYToutcont MPHYToutcont GPHYToutcont %OUTPUTS 
global FZOOoutcont  EZOOoutcont  MZOOoutcont 
global FDINoutcont  FPONoutcont
%...................................................................................

%%%%%%%%%%%%%%%%%
%STATE VARIABLES:
%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
% eps constraint added by Le Gland (21/11/2019)
phy  = max(V0(Iphy),eps);
zoo  = max(V0(Izoo),eps);
din  = max(V0(Idin),eps);
pon  = max(V0(Ipon),eps);
box  = max(V0(Ibox),eps);
%...................................................................................
%===================================================================================
%VERTICAL MIXING OF STATISTICAL MOMENTS:
%...................................................................................
if strcmp(keyPhysics,'not')
    %...............................................................................
    xave  = V0(Ixave) ;
    yave  = V0(Iyave) ;
    xxvar = V0(Ixxvar);
    yyvar = V0(Iyyvar);
    xycov = V0(Ixycov);
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    xave_star  = V0(Ixave) ;
    yave_star  = V0(Iyave) ;
    xxvar_star = V0(Ixxvar);
    yyvar_star = V0(Iyyvar);
    xycov_star = V0(Ixycov);
    %
    %...............................................................................
    xave  = (xave_star ./phy);
    yave  = (yave_star ./phy);
    xxvar = (xxvar_star./phy) - xave.^2   ;
    yyvar = (yyvar_star./phy) - yave.^2   ;
    xycov = (xycov_star./phy) - xave.*yave;
    %...............................................................................
end
% Protection against negative variances and out-of-range correlations
varmin = 10^(-6); % "eps" is not fit, because eps*eps = 0;
xxvar = max(varmin, xxvar);
yyvar = max(varmin, yyvar);

% Avoid correlation smaller than 1 and add correction term
dXYCOVdt_control = zeros(ndepths,1);
for i=1:ndepths
    %if sqrt(xxvar(i)*yyvar(i))-eps < abs(xycov(i))
    if sqrt(xxvar(i)*yyvar(i))*0.999 < abs(xycov(i))    % Stricter conditions to avoid bugs in the results (Le Gland, 20/07/2020)
        dXYCOVdt_control(i) = (sqrt(xxvar(i)*yyvar(i))*0.998 - abs(xycov(i)))*sign(xycov(i))/deltat;
        xycov(i) = 0.999*sqrt(xxvar(i)*yyvar(i)).*sign(xycov(i));
    end
end


%...................................................................................
%===================================================================================
%................................................................................... 
VNPZD0 = [phy;zoo;din;pon]; 
%................................................................................... 
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DAY OF SIMULATION ANT TIME COUNTER:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
[jday,newday] = SPEAD_1D_daycounter(jday,iTime); %
%...................................................................................
jcounter = floor(iTime/deltat); %For ode4.
% $$$ jcounter = floor(iTime/deltat) + 1; %For ode1.
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK FOR NEGATIVE CONCENTRATIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
Ineg = find(VNPZD0 < 0);
%........................................................................
if ~isempty(Ineg > 0)
    iTime
    disp(['P , Z , N , D']);
    wconcs = [phy,zoo,din,pon]
    wconcsNeg = VNPZD0(Ineg);
    disp('Error!!! there are NEGATIVE concentrations!')
    pause
end
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%
%MASS CONSERVATION:
%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
ntoti = sum([phy + zoo + din + pon + box]); %Checking mass conservation.
%...................................................................................
%%ydistmax = 1d-6; %original. 
ydistmax = 1d-3; %less strict level (for fast-numerical-solving)
ydist = abs(ntoti - ntot0);
%...................................................................................
%%if strcmp(keyNutrientSupply,'not') 
if abs(ydist) > ydistmax 
    masscheck_N = [iTime,jday,ntot0,ntoti,ydist]
    disp('Error!!! mass is NOT conserved!')
    pause
end 
%%end 
%...................................................................................
if mod(jday,10) == 0
    if strcmp(newday,'yes')
    masscheck = [iTime,jday,ntot0,ntoti,ydist]
    % disp('-------------------------------------------------')
    end
end
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MEAN TRAIT (j) OF THE GAUSSIAN DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
xj = xave; 
% xmj = xave; 
% sigmaxj = sqrt(xxvar);
%...................................................................................
yj = yave;
% ymj = yave;
% sigmayj = sqrt(yyvar);
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK FOR NEGATIVE VARIANCE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
Inegx = xxvar < 0;
Inegy = yyvar < 0;
%...................................................................................
xxvar(Inegx) = sqrt(eps);
yyvar(Inegy) = sqrt(eps);
%...................................................................................
Jnegx = find(xxvar < 0);
Jnegy = find(yyvar < 0);
%...................................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% $$$ if iTime > 9 
% $$$     x_ave_var_std2 = [xave,xvar,xstd.^2]
% $$$ end
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%...................................................................................
%if length(Jneg > 0)
if ~isempty(Jnegx > 0 | Jnegy > 0)
    iTime
    %x_ave_var_std2 = [xave,xvar,xstd.^2]
    xy_ave_var = [xave,yave,xxvar,yyvar,xycov];
    disp('Error!!! x_variance is negative!!!')
    pause 
end
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%
%GAUSSIAN DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%
tempj = temp(:,jday);

if sum(galfa ~= 1) ~= 0    
    fxyj         = zeros(ndepths,1)  ;
    sigma        = zeros(ndepths,2,2);
    invsigma     = zeros(ndepths,2,2);
    detsigma     = zeros(ndepths,1)  ;
    detsigmasqrt = zeros(ndepths,1)  ;
    % "vecj" is present in comments to show how some derivatives are computed.
    % However, it is equal to zero (all values and derivatives are taken around the mean value).
    %vecj         = zeros(ndepths,2)  ;
    for i=1:ndepths
       sigma(i,:,:)        = [xxvar(i), xycov(i); xycov(i), yyvar(i)];  % Variance matrix [trait^2]
       % This is only necessary in KTW mode (Le Gland, 26/11/2019)
       % invsigma(i,:,:)     = inv(squeeze(sigma(i,:,:)));                % Inverse of variance matrix [trait^(-2)]
       detsigma(i)         = (xxvar(i) .* yyvar(i) - xycov(i).^2);      % Determinant of variance matrix [trait^4]
       detsigmasqrt(i)     = sqrt(detsigma(i));                         % Square root of variance determinant [trait^2]
       %vecj(i,:)           = [xj(i)-xmj(i),yj(i)-ymj(i)];               % Trait vector compared to mean trait [trait]
       %vecjnow = squeeze(vecj(i,:));
    
       % fxyj(i) = (1.0 ./ (detsigmasqrt(i)*2*pi)) .* exp( -(1/2)*vecjnow*invsigma(i)*vecjnow' ); % 2D probability density function [trait^(-2)]
       % use \ instead of inverse matrix, to increase precision and decrease cost
       % fxyj(i) = (1.0 ./ (detsigmasqrt(i)*2*pi)) .* exp( -(1/2)*vecjnow*((squeeze(sigma(i,:,:)\vecjnow') );
       % Anyway, this should be zero
       fxyj(i) = (1.0 ./ (detsigmasqrt(i)*2*pi));
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PHYTOPLANKTON BIOMASS GAUSSIAN DISTRIBUTION:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ptot = phy;
    Pxyj = Ptot .* fxyj;
    Pxyjalfa = Pxyj.^galfa;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ANALYTICAL INTEGRATION OF PHYTOPLANKTON GASSIAN DISTRUBUTION:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===================================================================================
    %...................................................................................
    intPxyalfadxdy = (Ptot.^galfa ./ galfa) .* (detsigmasqrt*2*pi).^(1-galfa); % 2 traits 
    % intPxyalfadxdy = (Ptot.^galfa ./ sqrt(galfa)) .* (sigmaxj*sqrt(2*pi)).^(1-galfa); % KTW on size only
    
    %...................................................................................
    %===================================================================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %GRAZING FUNCTIONAL RESPONSE KTW:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===================================================================================
    %...................................................................................
    % Vmax = (gzmax*zoo);
    Vmax = (gzmax*zoo).*(Q10h.^((tempj - temp0)/10));
    Qswitchj = (Pxyjalfa ./ intPxyalfadxdy);
    Qfeeding = (Ptot.^gbeta ./ (Ptot.^gbeta + kgz^gbeta));
    %...................................................................................
    Gxyj = Qswitchj .* Qfeeding .* Vmax;
    %===================================================================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %BIOMASS SPECIFIC GRAZING FUNCTIONAL RESPONSE KTW:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %...................................................................................
    gxyj = Gxyj./Pxyj;
    %...................................................................................

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ANALYTICAL DERIVATIVES OF BIOMASS SPECIFIC GRAZING FUNCTIONAL RESPONSE KTW:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===================================================================================
    %...................................................................................
    d1gxydx   = zeros(ndepths,1);
    d1gxydy   = zeros(ndepths,1);
    d2gxydxdx = zeros(ndepths,1);
    d2gxydydy = zeros(ndepths,1);
    d2gxydxdy = zeros(ndepths,1);
    for i = 1:ndepths
       %d1gxydx(i)   = - gxyj(i) * (galfa(i) - 1) * (vecj(i,:) * squeeze(invsigma(i,1,:)) );
       %d1gxydy(i)   = - gxyj(i) * (galfa(i) - 1) * (vecj(i,:) * squeeze(invsigma(i,2,:)) );
       %d2gxydxdx(i) = gxyj(i) * (galfa(i) - 1) * ( (galfa(i) - 1) * ( vecj(i,:) * squeeze(invsigma(i,1,:)) )^2 - invsigma(i,1,1) );
       %d2gxydydy(i) = gxyj(i) * (galfa(i) - 1) * ( (galfa(i) - 1) * ( vecj(i,:) * squeeze(invsigma(i,2,:)) )^2 - invsigma(i,2,2) );
       %d2gxydxdy(i) = gxyj(i) * (galfa(i) - 1) * ( (galfa(i) - 1) * ( vecj(i,:) * squeeze(invsigma(i,1,:)) )*(vecj(i,:) * squeeze(invsigma(i,2,:)) ) - invsigma(i,1,2) );
       d2gxydxdx(i) = - gxyj(i) * (galfa(i) - 1) * invsigma(i,1,1);
       d2gxydydy(i) = - gxyj(i) * (galfa(i) - 1) * invsigma(i,2,2);
       d2gxydxdy(i) = - gxyj(i) * (galfa(i) - 1) * invsigma(i,1,2);
       % KTW on size only
       % d2gxydxdx(i) = gxyj(i) * (galfa(i) - 1) * ( (galfa(i) - 1) * ( vecj(i,:) * squeeze(invsigma(i,1,:)) )^2 - invsigma(i,1,1) );
       % d2gxydydy(i) = 0;
       % KTW on temperature only
       % d2gxydxdx(i) = 0;
       % d2gxydydy(i) = gxyj(i) * (galfa(i) - 1) * ( (galfa(i) - 1) * ( vecj(i,:) * squeeze(invsigma(i,2,:)) )^2 - invsigma(i,2,2) );
       % d2gxydxdy(i) = 0;
    end

% Simplify the case without KTW
elseif sum(galfa ~= 1) == 0
    Vmax = (gzmax*zoo).*(Q10h.^((tempj - temp0)/10));
    Qfeeding = (phy.^gbeta ./ (phy.^gbeta + kgz^gbeta));
    d1gxydx   = 0;
    d1gxydy   = 0;
    d2gxydxdx = 0;
    d2gxydydy = 0;
    d2gxydxdy = 0;
end
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%
%TURBULENT DIFFUSION:
%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
if ndepths > 1
kz  = KZ (:,jday);
end
%...................................................................................
%===================================================================================
%STATISTICAL MOMENTS:
if strcmp(keyPhysics,'yes')
%...................................................................................
DIFFxave_star  = zeros(ndepths,1);
DIFFxxvar_star = zeros(ndepths,1);
DIFFyave_star  = zeros(ndepths,1);
DIFFyyvar_star = zeros(ndepths,1);
DIFFxycov_star = zeros(ndepths,1);
%...................................................................................
if ndepths > 1
    [DIFFxave_star]  = SPEAD_1D_TurbulentDiffusion(xave_star,deltat,deltaz,kz,ndepths,'Implicit',100);
    [DIFFxxvar_star] = SPEAD_1D_TurbulentDiffusion(xxvar_star,deltat,deltaz,kz,ndepths,'Implicit',100);
    [DIFFyave_star]  = SPEAD_1D_TurbulentDiffusion(yave_star,deltat,deltaz,kz,ndepths,'Implicit',100);
    [DIFFyyvar_star] = SPEAD_1D_TurbulentDiffusion(yyvar_star,deltat,deltaz,kz,ndepths,'Implicit',100);
    [DIFFxycov_star] = SPEAD_1D_TurbulentDiffusion(xycov_star,deltat,deltaz,kz,ndepths,'Implicit',100);
end
%...................................................................................
end
%===================================================================================
%PLANKTON BIOMASSES:
%...................................................................................
DIFFphy = zeros(ndepths,1);
DIFFzoo = zeros(ndepths,1);
DIFFdin = zeros(ndepths,1);
DIFFpon = zeros(ndepths,1);
%...................................................................................
%pon
if ndepths > 1
    [DIFFphy] = SPEAD_1D_TurbulentDiffusion(phy,deltat,deltaz,kz,ndepths,'Implicit',100);
    [DIFFzoo] = SPEAD_1D_TurbulentDiffusion(zoo,deltat,deltaz,kz,ndepths,'Implicit',100);
    [DIFFdin] = SPEAD_1D_TurbulentDiffusion(din,deltat,deltaz,kz,ndepths,'Implicit',100);
    [DIFFpon] = SPEAD_1D_TurbulentDiffusion(pon,deltat,deltaz,kz,ndepths,'Implicit',100);
end
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%
%VERTICAL SINKING:
%%%%%%%%%%%%%%%%%%
%========================================================================
[ADVpon] = SPEAD_1D_SinkingAdvection(pon,deltaz,wsink,ndepths);
% Transform PON to DIN at bottom to avoid PON accumulation
ADVdin = zeros(ndepths,1);
ADVdin(end) = (wsink/deltaz)*pon(end);
%........................................................................
%========================================================================

%%%%%%%
%LIGHT:
%%%%%%%
%========================================================================
%........................................................................
jpar0 = parz0(jday); %Photo. Active. Radiation at the surface [W*m-2]
%........................................................................
jPAR  = jpar0*exp(-kw*zdepths(:)); %[W*m-2] PAR profile.
%........................................................................
%Qpar = (jPAR / Isat).*exp(1 - (jPAR/Isat)); %Phy light limitation [n.d.] values between 0 and 1.
% Use normalized Follows (2007) formula, with photoinhibition of large cells
Kpar   = log(InhFac+1) / Isat;
Kinhib = Kpar / InhFac; 
Fmax   = (Kpar+Kinhib)/Kpar * exp( -(Kinhib/Kpar) * log(Kinhib/(Kpar+Kinhib)) );  
% Qpar is normalized to have a maximum of 1 at Isat
Qpar = Fmax * (1 - exp(-Kpar*jPAR)) .* exp(-Kinhib*jPAR);
%........................................................................
%========================================================================
%------------------------------------------------------------------------
%NOTE: To keep the same DIN trade-off, both "mup" and "alp" ** must ** be 
%multiplied by the environmental limitation factor (eg. Qpar or Qsst) 
%Otherwise, if Qpar or Qsst only multiplies "mup", the optimal size for a 
%given DIN value will shift up and down instead of remaining always at 
%the same ESDphy value. 
%------------------------------------------------------------------------
%........................................................................
% $$$ alp = alp0*ones(ndepths,1); 
% $$$ mup = mup0*ones(ndepths,1);
%........................................................................
alp = alp0*Qpar; 
mup = mup0*Qpar;
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHYTOPLANKTON NUTRIENT UPTAKE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHYTOPLANKTON GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================================
%...............................................................................
knp = (mup./alp); 
%...............................................................................
%===============================================================================
%EXPONENTIAL UPTAKE AFFINITY: 
%...............................................................................
knxj = knp .* exp(aknp*xj); %Phy half-sat uptake as a function of cell size [mmolN*m-3]
%...............................................................................
%===============================================================================
%EXPONENTIAL GROWTH RATE AT REFERENCE TEMPERATURE:
%...............................................................................
muxj = mup .* exp(amup*xj);
%...............................................................................
%===============================================================================
%...............................................................................
lxj = (knxj ./ (knxj + din)); %[n.d.]
qxj = (din  ./ (din + knxj)); %[n.d.]
% Skewed response to temperature
q10 = Q10a.^((yj - temp0)/10);

if tempj < yj +5
    qyj = exp(0.2*(tempj - yj)) .* (yj + 5 - tempj)/5 .* q10;
else
    qyj = 0;
end
% qyj = exp(-(yj - tempj).^2 / (2*gammay^2)) .* q10;

%............................................................................... 
uxyj = muxj .* qxj .* qyj; %Uptake rate at mean temperature and mean size [d-1] 
%...............................................................................
%===============================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYTICAL DERIVATIVES OF GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================================
%MY DERIVATIONS FOR MICHAELIS MENTEN: 
%...............................................................................
d1qxdx = - qxj .* (aknp * lxj);
d1lxdx = - d1qxdx; 
%...............................................................................
d2qxdx = -aknp * (d1lxdx .* qxj + d1qxdx .* lxj); 
%...............................................................................
%===============================================================================
%FROM BINGZANG CHEN FOR MAXIMUM GROWTH RATE WITH SIZE: 
%...............................................................................
cff = amup;
%...............................................................................
d1muxdx = muxj.*cff;
d2muxdx = muxj.*cff.^2;
%...............................................................................
%===============================================================================
%...............................................................................
d1uxydx = qyj .* ( (d1muxdx .* qxj) + (muxj .* d1qxdx) );
%...............................................................................
d2uxydxdx = qyj .* ( ((d2muxdx .*    qxj) + (d1muxdx .* d1qxdx)) + ...
          ((d1muxdx .* d1qxdx) + (muxj    .* d2qxdx)) );
%...............................................................................
%===============================================================================
% New skewed response to temperature
if tempj < yj + 5
    d1uxydy = uxyj .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 );
    d2uxydydy = d1uxydy .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 ) - uxyj .* 1./(yj + 5 - tempj).^2;    
    d2uxydxdy = qyj .* ( (d1muxdx .* qxj) + (muxj .* d1qxdx) ) .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 );
else
    d1uxydy = 0;
    d2uxydydy = 0;
    d2uxydxdy = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ORDINARY DIFFERENTIAL EQUATIONS (ODEs):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
uxy = uxyj + (1/2)*(xxvar.*d2uxydxdx + yyvar.*d2uxydydy) + xycov.*d2uxydxdy;
%...................................................................................
gxy = (1./phy) .* Qfeeding .* Vmax;
%...................................................................................
%===================================================================================
%...................................................................................
Fphy = uxy .* phy;
Gphy = gxy .* phy;
Mphy = mp * phy .* (Q10h.^((tempj - temp0)/10)); % Temperature dependence on mortality
%...................................................................................
Fzoo = Gphy; %Zoo second production [mmolN * m-3 * d-1]
Ezoo = (1-betaz) * Fzoo; %Zoo exudation [mmolN * m-3 * d-1] 
Mzoo = mz * (zoo.^mpower) .* (Q10h.^((tempj - temp0)/10)); % Temperature dependence on mortality (Le Gland, 11/11/2019)
%...................................................................................
% Temperature-dependent, based on heterotrophic Q10 (Le Gland, 04/11/2019)
mdt = md*Q10h.^((tempj-temp0)/10);
Mpon = mdt.*pon;
%...................................................................................
%===================================================================================
%...................................................................................
dPHYdt =   Fphy - Gphy - Mphy; 
%...................................................................................
dZOOdt =   Fzoo - Ezoo - Mzoo; 
%...................................................................................
%dDINdt = - Fphy + epsPhy*Ephy + omePhy*Mphy + epsZoo*Ezoo + omeZoo*Mzoo + Mpon; 
dDINdt = - Fphy + omePhy*Mphy + epsZoo*Ezoo + omeZoo*Mzoo + Mpon;
%...................................................................................
%dPONdt = (1-epsPhy)*Ephy + (1-omePhy)*Mphy + (1-epsZoo)*Ezoo + (1-omeZoo)*Mzoo - Mpon;
dPONdt = (1-omePhy)*Mphy + (1-epsZoo)*Ezoo + (1-omeZoo)*Mzoo - Mpon;
%...................................................................................
%===================================================================================
%.......................................  ..........................................
dXAVEdt = xxvar .* (d1uxydx - d1gxydx) + xycov .* (d1uxydy - d1gxydy);
dYAVEdt = xycov .* (d1uxydx - d1gxydx) + yyvar .* (d1uxydy - d1gxydy);
%...................................................................................
dXXVARdt = (xxvar.^2) .* (d2uxydxdx - d2gxydxdx) + (xxvar.*xycov) .* (d2uxydxdy - d2gxydxdy) + ...
           (xycov.^2) .* (d2uxydydy - d2gxydydy) + numutx .* (2.*uxy); % (2.*uxyj + xxvar.*d2uxydxdx + ...
           %yyvar.*d2uxydydy + 2.*xycov.*d2uxydxdy);
dYYVARdt = (xycov.^2) .* (d2uxydxdx - d2gxydxdx) + (yyvar.*xycov) .* (d2uxydxdy - d2gxydxdy) + ...
           (yyvar.^2) .* (d2uxydydy - d2gxydydy) + numuty .* (2.*uxy); % (2.*uxyj + xxvar.*d2uxydxdx + ...
           %yyvar.*d2uxydydy + 2.*xycov.*d2uxydxdy);
dXYCOVdt = (xxvar.*xycov) .* (d2uxydxdx - d2gxydxdx) + (yyvar.*xycov) .* (d2uxydydy - d2gxydydy) + ...
           (xxvar.*yyvar + xycov.^2) .* (d2uxydxdy - d2gxydxdy) + dXYCOVdt_control;
%................................................................................... 
%===================================================================================
%...................................................................................
FPHYToutcont(:,jcounter) = Fphy;
GPHYToutcont(:,jcounter) = Gphy;
MPHYToutcont(:,jcounter) = Mphy;
%...................................................................................
FZOOoutcont(:,jcounter) = Fzoo;
EZOOoutcont(:,jcounter) = Ezoo;
MZOOoutcont(:,jcounter) = Mzoo;
%...................................................................................
FDINoutcont(:,jcounter) = omePhy*Mphy + epsZoo*Ezoo + omeZoo*Mzoo + Mpon; 
FPONoutcont(:,jcounter) = (1-omePhy)*Mphy + (1-epsZoo)*Ezoo + (1-omeZoo)*Mzoo;
%...................................................................................
%===================================================================================
%...................................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%DEBUGGING:
% $$$ %...................................................................................
% $$$ %%if mod(ti,10)==0 %show every 10 days
% $$$ masscheck_N = [iTime,jday,ntot0,ntoti,ydist]
% $$$ wGdx = [gx,gxj,d1gxdx,d2gxdx]
% $$$ wUdx = [ux,uxj,d1uxdx,d2uxdx,d3uxdx,d4uxdx]
% $$$ % $$$ wBIOs = [xave,xvar,xstd,phy,zoo,din,pon]
% $$$ % $$$ wODEs = [dXAVEdt,dXVARdt,dXSTDdt,dPHYdt,dZOOdt,dDINdt,dPONdt]
% $$$ wBIOs = [phy,zoo,din,pon]
% $$$ wODEs = [dPHYdt,dZOOdt,dDINdt,dPONdt]
% $$$ %...................................................................................
% $$$ dispPHYdot = [dPHYdt,Fphy,Ephy,Fzoo,Mphy]
% $$$ dispZOOdot = [dZOOdt,Fzoo,Ezoo,Mzoo]
% $$$ dispDINdot = [dDINdt,Fphy,epsPhy*Ephy,omePhy*Mphy,epsZoo*Ezoo,omeZoo*Mzoo,Mpon] 
% $$$ dispPONdot = [dPONdt,(1-epsPhy)*Ephy,(1-omePhy)*Mphy,(1-epsZoo)*Ezoo,(1-omeZoo)*Mzoo,Mpon] 
% $$$ %...................................................................................
% $$$ disp('-------------')
% $$$ pause
% $$$ %%end
% $$$ %...................................................................................
% $$$ % $$$ dispPHYdot = [dPHYdt,dPHYdtBis]
% $$$ % $$$ dispZOOdot = [dZOOdt,dZOOdtBis]
% $$$ % $$$ dispDINdot = [dDINdt,dDINdtBis]
% $$$ % $$$ dispPONdot = [dPONdt,dPONdtBis]
% $$$ % $$$ pause 
% $$$ %...................................................................................
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%===================================================================================

%................................................................................... 
if strcmp(keyPhysics,'yes')
dXAVE_STARdt = (dPHYdt .* xave) + (dXAVEdt .* phy);
dXXVAR_STARdt = (dPHYdt .* xxvar) + (dXXVARdt .* phy) + (dPHYdt .* xave .* xave) + 2*(dXAVEdt .* xave .* phy);
dYAVE_STARdt  =  (dPHYdt .* yave)  + (dYAVEdt .* phy);
dYYVAR_STARdt =  (dPHYdt .* yyvar) + (dYYVARdt .* phy) + (dPHYdt .* yave .* yave) + 2*(dYAVEdt .* yave .* phy);
dXYCOV_STARdt =  (dPHYdt .* xycov) + (dXYCOVdt .* phy) + (dPHYdt .* xave .* yave) + (dXAVEdt .* yave .* phy) + (dYAVEdt .* xave .* phy);
end

%...................................................................................
dBOXdt = dPHYdt + dZOOdt + dDINdt + dPONdt; %Virtual box to check mass conservation.
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%
%ADD PHYSICAL PROCESSES:
%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
phydot = dPHYdt + DIFFphy;
zoodot = dZOOdt + DIFFzoo;
dindot = dDINdt + DIFFdin + ADVdin; % Bottom remineralization (ADVdin)
pondot = dPONdt + DIFFpon + ADVpon; %Only PON has vertical sinking.
%...................................................................................
boxdot = dBOXdt; 
%...................................................................................
%===================================================================================
if strcmp(keyPhysics,'not')
    %...............................................................................
    xavedot  = dXAVEdt;
    xxvardot = dXXVARdt;
    yavedot  = dYAVEdt;
    yyvardot = dYYVARdt;
    xycovdot = dXYCOVdt;
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    xave_stardot  = dXAVE_STARdt  + DIFFxave_star;
    xxvar_stardot = dXXVAR_STARdt + DIFFxxvar_star;
    yave_stardot  = dYAVE_STARdt  + DIFFyave_star;
    yyvar_stardot = dYYVAR_STARdt + DIFFyyvar_star;
    xycov_stardot = dXYCOV_STARdt + DIFFxycov_star;
    %...............................................................................
end
%...................................................................................
%===================================================================================


%%%%%%%%%%
%STOCKAGE:
%%%%%%%%%%
% Done through globals. Unactivated by default. Perhaps there is a best way
% to do it.
%===================================================================================
%...................................................................................
todedotday(1,jday) = jday;
%...................................................................................
%Xavedotday(:,jday) = dXAVEdt;
%Xvardotday(:,jday) = dXVARdt;

%Xstddotday(:,jday) = dXSTDdt;
%...................................................................................
%d1UXdxday(:,jday) = d1uxdx;
%d1GXdxday(:,jday) = d1gxdx;
%...................................................................................
%d2UXdxday(:,jday) = d2uxdx;
%d2GXdxday(:,jday) = d2gxdx;
%...................................................................................
%UXday(:,jday) = ux;
%GXday(:,jday) = gx;
%...................................................................................
%===================================================================================
%...................................................................................
todedotout(1,jcounter) = jcounter*deltat;
%...................................................................................
%Xavedotout(:,jcounter) = dXAVEdt;
%Xvardotout(:,jcounter) = dXVARdt;
%Xstddotout(:,jcounter) = dXSTDdt;
%...................................................................................
%d1UXdxout(:,jcounter) = d1uxdx;
%d1GXdxout(:,jcounter) = d1gxdx;
%...................................................................................
%d2UXdxout(:,jcounter) = d2uxdx;
%d2GXdxout(:,jcounter) = d2gxdx;
%...................................................................................
%UXout(:,jcounter) = ux;
%GXout(:,jcounter) = gx;
%...................................................................................
%===================================================================================

% Case with 2 traits (Le Gland, 16/07/2019)
% Xavedotday(:,jday)  = dXAVEdt;
% Yavedotday(:,jday)  = dYAVEdt;
% XXvardotday(:,jday) = dXXVARdt;
% YYvardotday(:,jday) = dYYVARdt;
% XYcovdotday(:,jday) = dXYCOVdt;
% %...................................................................................
% UXYday(:,jday) = uxy;
% GXYday(:,jday) = gxy;
% %...................................................................................
% d1UXYdxday(:,jday) = d1uxydx;
% d1GXYdxday(:,jday) = d1gxydx;
% d1UXYdyday(:,jday) = d1uxydy;
% d1GXYdyday(:,jday) = d1gxydy;
% %...................................................................................
% d2UXYdxdxday(:,jday) = d2uxydxdx;
% d2GXYdxdxday(:,jday) = d2gxydxdx;
% d2UXYdydyday(:,jday) = d2uxydydy;
% d2GXYdydyday(:,jday) = d2gxydydy;
% d2UXYdxdyday(:,jday) = d2uxydxdy;
% d2GXYdxdyday(:,jday) = d2gxydxdy;
% %...................................................................................
% %===================================================================================
% %...................................................................................
% Xavedotout(:,jcounter)  = dXAVEdt;
% Yavedotout(:,jcounter)  = dYAVEdt;
% XXvardotout(:,jcounter) = dXXVARdt;
% YYvardotout(:,jcounter) = dYYVARdt;
% XYcovdotout(:,jcounter) = dXYCOVdt;
% %...................................................................................
UXYout(:,jcounter) = uxy;
GXYout(:,jcounter) = gxy;
% %...................................................................................
% d1UXYdxout(:,jcounter) = d1uxydx;
% d1GXYdxout(:,jcounter) = d1gxydx;
% d1UXYdyout(:,jcounter) = d1uxydy;
% d1GXYdyout(:,jcounter) = d1gxydy;
% %...................................................................................
% d2UXYdxdxout(:,jcounter) = d2uxydxdx;
% d2GXYdxdxout(:,jcounter) = d2gxydxdx;
% d2UXYdydyout(:,jcounter) = d2uxydydy;
% d2GXYdydyout(:,jcounter) = d2gxydydy;
% d2UXYdxdyout(:,jcounter) = d2uxydxdy;
% d2GXYdxdyout(:,jcounter) = d2gxydxdy;
%...................................................................................
%===================================================================================
%...................................................................................

%%%%%%%%
%OUTPUT:
%%%%%%%%
%===================================================================================
%...................................................................................
if strcmp(keyPhysics,'not')
    %...............................................................................
    Vdot = [phydot;zoodot;dindot;pondot;boxdot;xavedot;yavedot;xxvardot;yyvardot;xycovdot];
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    Vdot = [phydot;zoodot;dindot;pondot;boxdot;xave_stardot;yave_stardot;xxvar_stardot;yyvar_stardot;xycov_stardot];
    %...............................................................................
end

return
%...................................................................................
%===================================================================================
%***********************************************************************************
