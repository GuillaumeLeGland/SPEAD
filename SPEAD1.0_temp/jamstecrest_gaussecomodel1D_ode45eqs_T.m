function [Vdot] = jamstecrest_gaussecomodel1D_ode45eqs_T(iTime,V0)
global galfa gbeta
% global gzmax kgz mz betaz betap mpower 
global gzmax kgz mz betaz betap_max mpower %betap_xmax betap_xrange  
global mp Isat InhFac numuty % numutx % Replace numut by 2 different mutation rates (Le Gland, 10/09/2019) 
global alp0 mup0 % knp0 
%global amup aknp bmup % aalp 
%global Q10 
global Q10a Q10h % Distinct partition coefficients for auto and heterotrophic processes (Le Gland, 31/10/2019)
global ntot0 
global temp0 temp % Use temperature at all depths (Le Gland, 15/10/2019) 
global jcounter 
%global iTimeode 
global  keyPhysics %keyAssimConstant %keyTraitAxis keySinking % keyAssimConstant added by Le Gland, 03/10/2019
global Sdin Drate % Mdin 
global epsPhy omePhy epsZoo omeZoo md 
global deltat %t0 ndays nyear tmax tspan 
global zdepths ndepths deltaz
global jday % jjday
%global Iave Ivar Istd %continuous model
global Iyave_T Iyyvar_T
global Iphy Izoo Idin Ipon Ibox %continuous model
global KZ % KZI
global parz0
global kw wsink % kp
% global keyNutrientSupply 
%...................................................................................
% global Xavedotday Xvardotday Xstddotday %OUTPUTS 
global UYday GYday 
% global d1UXdxday d1GXdxday 
% global d2UXdxday d2GXdxday 
global todedotday 
%...................................................................................
% global Xavedotout Xvardotout Xstddotout %OUTPUTS 
% global UXout GXout 
% global d1UXdxout d1GXdxout 
% global d2UXdxout d2GXdxout 
global todedotout 
% I do not know if these outputs are really useful (Le Gland, 25/11/2019)
% Outputs in the case with 2 traits (Le Gland, 16/07/2019)
% global Xavedotday Yavedotday XXvardotday YYvardotday XYcovdotday %OUTPUTS 
% global UXYday GXYday 
% global d1UXYdxday d1UXYdyday d1GXYdxday d1GXYdyday 
% global d2UXYdxdxday d2UXYdydyday d2UXYdxdyday d2GXYdxdxday d2GXYdydyday d2GXYdxdyday
% Outputs in the case with 2 traits (Le Gland, 16/07/2019)
% global Xavedotout Yavedotout XXvardotout YYvardotout XYcovdotout %OUTPUTS 
% global UXYout GXYout 
% global d1UXYdxout d1UXYdyout d1GXYdxout d1GXYdyout 
% global d2UXYdxdxout d2UXYdydyout d2UXYdxdyout d2GXYdxdxout d2GXYdydyout d2GXYdxdyout 
%...................................................................................
global FPHYToutcont EPHYToutcont MPHYToutcont GPHYToutcont %OUTPUTS 
global FZOOoutcont  EZOOoutcont  MZOOoutcont 
global FDINoutcont  FPONoutcont
%...................................................................................
% global DIFFxaveout DIFFxvarout DIFFxstdout 
%...................................................................................

%%%%%%%%%%%%%%%%%
%STATE VARIABLES:
%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
% eps constraint added by Le Gland (21/11/2019)
phy  = max(V0(Iphy),sqrt(eps));
zoo  = max(V0(Izoo),sqrt(eps));
din  = max(V0(Idin),sqrt(eps));
pon  = max(V0(Ipon),sqrt(eps));
box  = max(V0(Ibox),sqrt(eps));
%...................................................................................
%===================================================================================
%VERTICAL MIXING OF STATISTICAL MOMENTS:
%...................................................................................
if strcmp(keyPhysics,'not')
    %...............................................................................
    yave  = V0(Iyave_T) ;
    yyvar = V0(Iyyvar_T);
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    yave_star  = V0(Iyave_T) ;
    yyvar_star = V0(Iyyvar_T);
    %...............................................................................
    yave  = (yave_star ./phy);
    yyvar = (yyvar_star./phy) - yave.^2;
    %...............................................................................
end
% Protection against negative variances and out-of-range correlations (Le Gland, 21/11/2019)
yyvar = max(10*sqrt(eps), yyvar); 

%...................................................................................
%===================================================================================
%................................................................................... 
VNPZD0 = [phy;zoo;din;pon]; 
%...................................................................................
%VSTAT0 = [xave,xxvar];
%................................................................................... 
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DAY OF SIMULATION ANT TIME COUNTER:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
[jday,newday] = jamstecrest_daycounter(jday,iTime); % I simplify the function (Le Gland, 13/09/2019)
%...................................................................................
%jcounter1 = jcounter; %Previou-s jcounter.
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
    wdilut = [Drate(jcounter),Sdin(jcounter)]
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
ntoti = sum(phy + zoo + din + pon + box); %Checking mass conservation.
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
if mod(jday,10) == 0 % jjday is useless and can be replaced by jday (Le Gland, 13/09/2019)
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
yj = yave; 
%ymj = yave; 
%...................................................................................
sigmayj = sqrt(yyvar);
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK FOR NEGATIVE VARIANCE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
Inegx = yyvar < 0;
%...................................................................................
yyvar(Inegx) = sqrt(eps);
%...................................................................................
Jnegx = find(yyvar < 0);
%...................................................................................
if ~isempty(Jnegx > 0)
    iTime
    xy_ave_var = [xave,xxvar]
    disp('Error!!! x_variance is negative!!!')
    pause 
end
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%
%GAUSSIAN DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%
%fxj = (1.0 ./ (sigmaxj * sqrt(2*pi))) .* exp( -(xj - xmj).^2 ./ (2*sigmaxj.^2) ); %Okay.
fyj = (1.0 ./ (sigmayj * sqrt(2*pi))); % since xj=xmj (Le Gland, 22/11/2019) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHYTOPLANKTON BIOMASS GAUSSIAN DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ptot = phy;
Pyj = Ptot .* fyj;
Pyjalfa = Pyj.^galfa;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYTICAL INTEGRATION OF PHYTOPLANKTON GASSIAN DISTRUBUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%intPxdx = Ptot; % Useless, since Ptot can be used instead (Le Gland, 25/11/2019)
%...................................................................................
intPyalfady = (Ptot.^galfa ./ sqrt(galfa)) .* (sigmayj*sqrt(2*pi)).^(1-galfa);
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYTICAL DERIVATIVES OF PHYTOPLANKTON GAUSSIAN DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%dPxdx = -(xj - xmj) .* (Pxj ./ sigmaxj.^2);
%dPxdx = 0; % Value at xj=xmj (Le Gland, 25/11/2019)
%...................................................................................
%d2Pxdx = Pxj .* ( ((xj - xmj).^2 ./ sigmaxj.^4) - (1.0 ./ sigmaxj.^2) );
%d2Pxdx = -Pxj ./ (sigmaxj.^2); % Value at xj=xmj (Le Gland, 25/11/2019)
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYTICAL DERIVATIVES OF PHYTOPLANKTON-ALFA GAUSSIAN DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%dPxalfadx = galfa .* (Pxj.^galfa) .* (- (xj - xmj) ./ sigmaxj.^2); %Okay!!!
%dPxalfadx = 0; % Since xj=xmj (Le Gland, 22/11/2019)
%...................................................................................
%d2Pxalfadx = galfa.*Pxjalfa .* (galfa.*((xj - xmj).^2 ./ sigmaxj.^4) - (1.0./sigmaxj.^2));
%d2Pxalfadx = -galfa.*Pxjalfa .* (1.0./sigmaxj.^2);
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GRAZING FUNCTIONAL RESPONSE KTW:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
tempj = temp(:,jday); % (Le Gland, 04/11/2019)
Vmax = (gzmax*zoo).*(Q10h.^((tempj - temp0)/10));
Qswitchj = (Pyjalfa ./ intPyalfady);
Qfeeding = (Ptot.^gbeta ./ (Ptot.^gbeta + kgz^gbeta));
%...................................................................................
Gyj = Qswitchj .* Qfeeding .* Vmax;
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYTICAL DERIVATIVES OF GRAZING KTW FUNCTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
% This part is useless, and not used later (Le Gland, 25/11/2019)
% A = (sigmaxj*sqrt(2*pi)).^(1-galfa) .* (1./sqrt(galfa));
% size only
%...................................................................................
% G = (Qfeeding .* Vmax); 
%...................................................................................
% B = (G./A); 
%...................................................................................
% d1Gxdx = galfa .* fxj.^galfa .* (- (xj - xmj) ./ sigmaxj.^2) .* B;
%...................................................................................
% d2Gxdx = galfa .* fxj.^galfa .* (galfa .* ((xj - xmj).^2 ./ sigmaxj.^4) - (1./sigmaxj.^2)) .* B;
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BIOMASS SPECIFIC GRAZING FUNCTIONAL RESPONSE KTW:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...................................................................................
gyj = Gyj./Pyj; %[d-1] 
%...................................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYTICAL DERIVATIVES OF BIOMASS SPECIFIC GRAZING FUNCTIONAL RESPONSE KTW:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
% d1gxdx = gxj .* (galfa - 1) .* (-(xj - xmj) ./ sigmaxj.^2 ); 
d1gydy = 0; % Value at xj=xmj (Le Gland, 25/11/2019)
%...................................................................................
% d2gxdx = gxj .* (galfa - 1) .* ( (galfa - 1) .* (-(xj - xmj) ./ sigmaxj.^2).^2 - (1./sigmaxj.^2) );
d2gydydy = -gyj .* (galfa-1) ./ (sigmayj.^2); % Value at xj=xmj (Le Gland, 25/11/2019)
%...................................................................................
%===================================================================================
    
% else
%     fxj = (xj == xmj);
%     Vmax = (gzmax*zoo);
%     Qfeeding = ( phy.^gbeta ./ ( phy.^gbeta + kgz^gbeta ) );
%     Gxj = Qfeeding .* Vmax;
%     d1gxdx = 0;
%     d2gxdx = 0;
% end

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
    DIFFyave_star  = zeros(ndepths,1);
    DIFFyyvar_star = zeros(ndepths,1);
    %...................................................................................
    if ndepths > 1
        [DIFFyave_star]  = jamstecrest_TurbulentDiffusion(yave_star,deltat,deltaz,kz,ndepths,'Reflectante');
        [DIFFyyvar_star] = jamstecrest_TurbulentDiffusion(yyvar_star,deltat,deltaz,kz,ndepths,'Reflectante');
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
    [DIFFphy] = jamstecrest_TurbulentDiffusion(phy,deltat,deltaz,kz,ndepths,'Reflectante');
    [DIFFzoo] = jamstecrest_TurbulentDiffusion(zoo,deltat,deltaz,kz,ndepths,'Reflectante');
    [DIFFdin] = jamstecrest_TurbulentDiffusion(din,deltat,deltaz,kz,ndepths,'Reflectante');
    [DIFFpon] = jamstecrest_TurbulentDiffusion(pon,deltat,deltaz,kz,ndepths,'Reflectante');
end
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%
%VERTICAL SINKING:
%%%%%%%%%%%%%%%%%%
%========================================================================
% if strcmp(keyPhysics,'not')
%     %....................................................................
%     ADVxave  = zeros(ndepths,1);
%     ADVxxvar = zeros(ndepths,1);
%     %....................................................................
% elseif strcmp(keyPhysics,'yes')
%     %....................................................................
%     ADVxave_star  = zeros(ndepths,1);
%     ADVxxvar_star = zeros(ndepths,1);
%     %....................................................................
% end
%........................................................................
%========================================================================
%........................................................................
%ADVphy = zeros(ndepths,1);
%ADVzoo = zeros(ndepths,1);
%ADVdin = zeros(ndepths,1);
%ADVpon = zeros(ndepths,1);
%........................................................................
[ADVpon] = jamstecrest_SinkingAdvection(pon,deltaz,wsink,ndepths);
% Transform PON to DIN at bottom to avoid PON accumulation (Le Gland, 11/12/2019)
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
% Qpar = (jPAR / Isat).*exp(1 - (jPAR/Isat)); %Phy light limitation [n.d.] values between 0 and 1.
% Suppress photoinhibition and uses a smaller and more realistic Isat (Le
% Gland, 31/10/2019)
% Qpar = tanh(jPAR/Isat);
% Use normalized Follows (2007) formula, with photoinhibition of large cells (Le Gland, 13/01/2020)
Kpar   = log(InhFac+1) / Isat;
Kinhib = Kpar / InhFac; 
Fmax   = (Kpar+Kinhib)/Kpar * exp( -(Kinhib/Kpar) * log(Kinhib/(Kpar+Kinhib)) );  
% Qpar is normalized to have a maximum of 1 at Isat (Le Gland, 13/01/2020)
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
%knxj = knp .* exp(aknp*xj); %Phy half-sat uptake as a function of cell size [mmolN*m-3]
%knxj = exp(xj); % Case where trait is half-saturation (Le Gland, 26/06/2019)
%kn = knp * exp(-0.5);% Same constant kn for all species (Le Gland, 25/11/2019)
kn = knp * exp(0); % Intent to better adapt winter and depth (Le Gland, 28/11/219)
%...............................................................................
%===============================================================================
%EXPONENTIAL MAXIMUM GROWTH RATE: 
%...............................................................................
%muxj = mup .* exp(amup*xj); %Phy maximum grazing rate as a function of cell size [d-1]
%mu = mup * exp(-0.25);
mu = mup * exp(0);
%...............................................................................
%===============================================================================
%UNIMODAL MAXIMUM GROWTH RATE AT REFERENCE TEMPERATURE (15 DEGREES):
%...............................................................................
%muxj = mup .* exp(amup*xj + bmup*xj.^2);
%muxj = mup .* exp(amup*(xj-log(knp))/aknp); % Case where trait is half-saturation (Le Gland, 26/06/2019)
%...............................................................................
%===============================================================================
%...............................................................................
%lxj = (knxj ./ (knxj + din)); %[n.d.]
%qx = (din  ./ (din + kn)); %[n.d.]

% New approach (20/07/2020): kn is not fixed for all phenotypes but always equal to 
% the local DIN concentration (less arbitrary, but also less similar to other models)
qx = 1/2; % din / (din + din)
mu = mup .* sqrt(din);
%...............................................................................

% New skewed response to temperature (Le Gland, 23/10/2019)
% q10 = Q10.^((yj - sst0)/10);
q10 = Q10a.^((yj - temp0)/10); % (Le Gland, 04/11/2019)
% q10 = Q10a.^(tempj-temp0); % q10 for single trait (Le Gland, 25/11/2019)

if tempj < yj +5
    qyj = exp(0.2*(tempj - yj)) .* (yj + 5 - tempj)/5 .* q10;
else
    qyj = 0;
end
%qyj = exp(-(yj - tempj).^2 / (2*gammay^2)) .* q10;

%............................................................................... 
%uxyj = muxj .* qxj .* qyj; %Uptake rate at mean temperature and mean size [d-1] 
uyj = mu .* qx .* qyj; % Uptke rate at mean size (Le Gland, 25/11/2019)
%...............................................................................
%===============================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYTICAL DERIVATIVES OF GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================================
%MY DERIVATIONS FOR MICHAELIS MENTEN: 
%...............................................................................
%d1qxdx = - qxj .* (aknp * lxj);
%d1qxdx = - qxj .* lxj; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%d1lxdx = - d1qxdx; 
%...............................................................................
%d2qxdx = -aknp * (d1lxdx .* qxj + d1qxdx .* lxj); 
%d2qxdx = d1lxdx .* qxj + d1qxdx .* lxj; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%d2lxdx = - d2qxdx; 
%...............................................................................
% Order 3 and 4 are useless (Le Gland, 25/11/2019)
%d3qxdx = -aknp * (d2lxdx .* qxj + d2qxdx .* lxj + 2 * d1lxdx .* d1qxdx); 
%d3qxdx = d2lxdx .* qxj + d2qxdx .* lxj + 2 * d1lxdx .* d1qxdx; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%d3lxdx = - d3qxdx; 
%...............................................................................
%d4qxdx = -aknp * ((d3lxdx .*   qxj  + d2lxdx .* d1qxdx) + ...
%           (d3qxdx .*   lxj  + d2qxdx .* d1lxdx) + ...
%           (d2lxdx .* d1qxdx + d1lxdx .* d2qxdx) * 2);
%d4qxdx = ((d3lxdx .*   qxj  + d2lxdx .* d1qxdx) + ...
%              (d3qxdx .*   lxj  + d2qxdx .* d1lxdx) + ...
%              (d2lxdx .* d1qxdx + d1lxdx .* d2qxdx) * 2); % Case where trait is half-saturation (Le Gland, 26/06/2019)
              
%...............................................................................
%===============================================================================
%FROM BINGZANG CHEN FOR MICHAELIS MENTEN: 
%...............................................................................
% $$$     d1qxdx_bis = -aknp*knxj.*qxj.^2;
% $$$     d2qxdx_bis = -aknp^2*din.*knxj.*(1./(din+knxj).^2 - 2*knxj./(din+knxj).^3);
% $$$     d3qxdx_bis =  aknp^3*din.*knxj.*(2*din.*knxj-(knxj-din).^2)./(knxj+din).^4;
% $$$     d4qxdx_bis =  aknp^4*din.*knxj.*(11*knxj.*din.*(din-knxj)+knxj.^3-din.^3)./(din+knxj).^5; 
    %...............................................................................
% $$$     d1lxdx = - d1qxdx_bis; 
% $$$     d2lxdx = - d2qxdx_bis; 
% $$$     d3lxdx = - d3qxdx_bis; 
%...............................................................................
%===============================================================================
%MY DERIVATIONS FOR NUTRIENT UPTAKE USING EXPONENTIAL GROWTH RATE WITH SIZE: 
%...............................................................................
% $$$     gff = (amup - aknp * lxj);
% $$$     %...............................................................................
% $$$     d1uxdx =   uxj  .*  gff; 
% $$$     %...............................................................................
% $$$     d2uxdx =   uxj  .* (gff.^2 - aknp^2 * lxj .* qxj); 
% $$$     %...............................................................................
% $$$     d3uxdx = d2uxdx .*  gff - 2*aknp * (d1uxdx .* d1lxdx) - aknp * (uxj .* d2lxdx); 
% $$$     %...............................................................................
% $$$     d4uxdx = d3uxdx .*  gff - 3*aknp * (d2uxdx .* d1lxdx + d1uxdx .* d2lxdx) - aknp * (uxj .* d3lxdx); 
%...............................................................................
%===============================================================================
%FROM BINGZANG CHEN FOR UNIMODAL MAXIMUM GROWTH RATE WITH SIZE: 
%...............................................................................
%cff = (amup + 2*bmup*xj); 
% cff = amup / aknp; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%...............................................................................
%d1muxdx = muxj.*cff;
%d2muxdx = muxj*2*bmup + muxj.*cff.^2;
% Order 3 and 4 are useless (Le Gland, 25/11/2019)
%d3muxdx = (2*bmup+cff.^2).*d1muxdx + 4.*bmup*muxj.*cff;
%d4muxdx = d1muxdx*8*bmup.*cff + (2*bmup+cff.^2).*d2muxdx + 8*bmup.^2*muxj;
%d2muxdx = muxj.*cff.^2; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%d3muxdx = muxj.*cff.^3; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%d4muxdx = muxj.*cff.^4; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%...............................................................................
%===============================================================================
%FROM BINGZANG CHEN FOR NUTRIENT UPTAKE USING UNIMODAL MAXIMUM GROWTH RATE: 
%...............................................................................
%d1uxdx = qyj .* ( (d1muxdx .* qxj) + (muxj .* d1qxdx) );
%...............................................................................
%d2uxdxdx = qyj .* ( ((d2muxdx .*    qxj) + (d1muxdx .* d1qxdx)) + ...
%          ((d1muxdx .* d1qxdx) + (muxj    .* d2qxdx)) );
%...............................................................................
% d3uxydx = qyj .* ( ((d3muxdx .*    qxj) + (d2muxdx .* d1qxdx)) + ...
%           ((d2muxdx .* d1qxdx) + (d1muxdx .* d2qxdx)) + ...
%           ((d2muxdx .* d1qxdx) + (d1muxdx .* d2qxdx)) + ...
%           ((d1muxdx .* d2qxdx) + (muxj    .* d3qxdx)) );
%...............................................................................
% d4uxydx = qyj .* ( ((d4muxdx .*    qxj) + (d3muxdx .* d1qxdx)) + ...
%           ((d3muxdx .* d1qxdx) + (d2muxdx .* d2qxdx)) + ...
%           ((d3muxdx .* d1qxdx) + (d2muxdx .* d2qxdx)) + ...
%           ((d2muxdx .* d2qxdx) + (d1muxdx .* d3qxdx)) + ...
%           ((d3muxdx .* d1qxdx) + (d2muxdx .* d2qxdx)) + ...
%           ((d2muxdx .* d2qxdx) + (d1muxdx .* d3qxdx)) + ...
%           ((d2muxdx .* d2qxdx) + (d1muxdx .* d3qxdx)) + ...
%           ((d1muxdx .* d3qxdx) + (muxj    .* d4qxdx)) );
%...............................................................................
%===============================================================================

% d1uydy = uyj .* (  -(yj - tempj)/gammay^2);
% d2uydydy = uyj .* ( (-(yj - tempj)/gammay^2).^2 - (1/gammay^2) );
% Third and fourth derivatives added, yj transformed to xj (Le Gland, 25/04/2019)
% d3uxydy = uxyj .* ( (-(yj - sstj)/gammay^2).^3 + (yj - sstj)/gammay^4 + 2*(yj - sstj)/gammay^4 );
% d4uxydy = uxyj .* ( (-(yj - sstj)/gammay^2).^4 - 3*((yj - sstj).^2)/gammay^6 - 3*((yj - sstj).^2)/gammay^6 + (3/gammay^4) );

% Cross derivative d2uxydxdy is necessary for further computations
% (Le Gland, 12/07/2019)
% d2uxydxdy = qyj .* ( (d1muxdx .* qxj) + (muxj .* d1qxdx) ) .* (  -(yj - tempj)/gammay^2);

% New symmetrical response (Le Gland, 15/11/2019)
% d1uxydy = uxyj .* (  -(yj - tempj)/gammay^2 + log(Q10a)/10 );
% d2uxydydy = uxyj .* ( (-(yj - tempj)/gammay^2 + log(Q10a)/10 ).^2 - (1/gammay^2) );
% d2uxydxdy = qyj .* ( (d1muxdx .* qxj) + (muxj .* d1qxdx) ) .* (  -(yj - tempj)/gammay^2 + log(Q10a)/10 );

% New skewed response to temperature (Le Gland, 23/10/2019)
if tempj < yj + 5
%     d1uxydy = uxyj .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 );
%     d2uxydydy = d1uxydy .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 ) - uxyj .* 1./(yj + 5 - tempj).^2; 
    d1uydy = uyj .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 );
    d2uydydy = d1uydy .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 ) - uyj .* 1./(yj + 5 - tempj).^2;
%     d2uxydxdy = qyj .* ( (d1muxdx .* qxj) + (muxj .* d1qxdx) ) .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 );
else
    d1uydy = 0;
    d2uydydy = 0;
%     d2uxydxdy = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ORDINARY DIFFERENTIAL EQUATIONS (ODEs):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
uy = uyj + (1/2)*(sigmayj.^2).*d2uydydy; %Community upake rate gain.

%Take exudation into account (Le Gland, 03/10/2019)
betap = betap_max;
ey = (1-betap)*uy;
d1eydy   = (1-betap).*d1uydy;
d2eydydy = (1-betap).*d2uydydy;

%...................................................................................
%===================================================================================
%...................................................................................
% $$$ gx = gxj + (1/2)*(sigmaxj.^2).*(d2gxdx); %Community grazing rate loss [d-1]. %WRONG!!!!
%...................................................................................
gy = (1./phy) .* Qfeeding .* Vmax; %Community grazing rate loss [d-1] %OKAY (NOTE: Ptot = phy)
% Test with grazing equal to growth (Le Gland, 09/05/2019)
% gx = uxj + (1/2)*(sigmaxj.^2).*d2uxdx;
%...................................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% $$$ disp('*** pto3 ***')
% $$$ whos phy Qfeeding Vmax 
% $$$ gx001 = gxj 
% $$$ gx002 = (1./phy) .* Qfeeding .* Vmax 
% $$$ pause 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%===================================================================================
%...................................................................................
Fphy = uy .* phy; %Phy primary production [mmolN * m-3 * d-1]
Gphy = gy .* phy; %Phy grazing mortality [mmolN * m-3 * d-1]
% Mphy = mp  * phy; %Phy natural mortality [mmolN * m-3 * d-1] %Linear.
Mphy = mp * phy .* (Q10h.^((tempj - temp0)/10)); % Temperature dependence on mortality (Le Gland, 11/11/2019)
% Mphy = mp  * (phy.*phy); %Phy natural mortality [mmolN * m-3 * d-1] %Quadratic.
% Mphy = mp * (phy.*phy) .* (Q10h.^((tempj - temp0)/10)); % Quadratic
% Test with grazing equal to growth (Le Gland, 09/05/2019)
% Mphy = 0;
%...................................................................................
Fzoo = Gphy; %Zoo second production [mmolN * m-3 * d-1]
Ezoo = (1-betaz) * Fzoo; %Zoo exudation [mmolN * m-3 * d-1] 
%Ephy = (1-betap) * Fphy; %Phy exudation [mmolN * m-3 * d-1]
Ephy = ey .* phy; % Variable exudation rate (Le Gland, 03/10/2019)
% Test with grazing equal to growth (Le Gland, 09/05/2019)
% Ephy = 0;
% Mzoo = mz * (zoo.^mpower); %Zoo natural mortality [mmolN * m-3 * d-1]
Mzoo = mz * (zoo.^mpower) .* (Q10h.^((tempj - temp0)/10)); % Temperature dependence on mortality (Le Gland, 11/11/2019)
%...................................................................................
%Mpon = md*pon; %Det degradation to nutrients [mmolN * m-3 * d-1]
% Temperature-dependent, based on heterotrophic Q10 (Le Gland, 04/11/2019)
mdt = md*Q10h.^((tempj-temp0)/10);
Mpon = mdt.*pon;
%...................................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% $$$ disp('*** pto2 ***')
% $$$ % $$$ gxj,
% $$$ % $$$ (1/2)*(sigmaxj.^2).*(d2gxdx)
% $$$ gx 
% $$$ phy 
% $$$ Fzoo 
% $$$ pause
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%===================================================================================
%-----------------------------------------------------------------------------------
% dsigma^2/dt = 2*sigma * dsigma/dt --> dsigma/dt = (1/(2*sigma)) * dsigma^2/dt
%-----------------------------------------------------------------------------------
%===================================================================================
%...................................................................................
dPHYdt =   Fphy - Ephy - Gphy - Mphy; 
%...................................................................................
dZOOdt =   Fzoo - Ezoo - Mzoo; 
%...................................................................................
dDINdt = - Fphy + epsPhy*Ephy + omePhy*Mphy + epsZoo*Ezoo + omeZoo*Mzoo + Mpon; 
%...................................................................................
dPONdt = (1-epsPhy)*Ephy + (1-omePhy)*Mphy + (1-epsZoo)*Ezoo + (1-omeZoo)*Mzoo - Mpon;
%...................................................................................
%===================================================================================
%...................................................................................
dYAVEdt = yyvar .* (d1uydy - d1gydy - d1eydy);
%....................................................................................
dYYVARdt = (yyvar.^2) .* (d2uydydy - d2gydydy - d2eydydy) + numuty .* (2.*(uy-ey)); %
%................................................................................... 
%===================================================================================
%...................................................................................
FPHYToutcont(:,jcounter) = Fphy;
EPHYToutcont(:,jcounter) = Ephy;
GPHYToutcont(:,jcounter) = Gphy;
MPHYToutcont(:,jcounter) = Mphy;
%...................................................................................
FZOOoutcont(:,jcounter) = Fzoo;
EZOOoutcont(:,jcounter) = Ezoo;
MZOOoutcont(:,jcounter) = Mzoo;
%...................................................................................
FDINoutcont(:,jcounter) = epsPhy*Ephy + omePhy*Mphy + epsZoo*Ezoo + omeZoo*Mzoo + Mpon; 
FPONoutcont(:,jcounter) = (1-epsPhy)*Ephy + (1-omePhy)*Mphy + (1-epsZoo)*Ezoo + (1-omeZoo)*Mzoo;
%...................................................................................
%===================================================================================
%...................................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADD EXTERNAL NUTRIENT SOURCES AND ALL STATE VARIABLES SINKS BY DILUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
Dphy = Drate(jcounter)*phy; %Dilution of phytoplankton [mmolN*m-3*d-1]
Dzoo = Drate(jcounter)*zoo;
Ddin = Drate(jcounter)*din;
Dpon = Drate(jcounter)*pon;
%...................................................................................
Sno3 = Sdin(:,jcounter); %Supply of nutrients into the system [mmolN*m-3*d-1]
%...................................................................................
%%Dzoo = zeros(ndepths,1); %Remove washout dilution of zooplanton.
%...................................................................................
%===================================================================================
%...................................................................................

dPHYdt = dPHYdt - Dphy;
dZOOdt = dZOOdt - Dzoo;
dDINdt = dDINdt - Ddin + Sno3; %Only DIN has an input flux source.
dPONdt = dPONdt - Dpon;

%................................................................................... 
if strcmp(keyPhysics,'yes')
    dYAVE_STARdt = (dPHYdt .* yave) + (dYAVEdt .* phy);
    % dXVAR_STARdt = (dPHYdt .* xvar) + (dXVARdt .* phy);
    % Change by Le Gland (23/04/2019)
    dYYVAR_STARdt = (dPHYdt .* yyvar) + (dYYVARdt .* phy) + (dPHYdt .* yave .* yave) + 2*(dYAVEdt .* yave .* phy);
end

%...................................................................................
dBOXdt = -Sno3 + Dphy + Dzoo + Ddin + Dpon; %Virtual box to check mass conservation.
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%
%ADD PHYSICAL PROCESSES:
%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
phydot = dPHYdt + DIFFphy;
zoodot = dZOOdt + DIFFzoo;
dindot = dDINdt + DIFFdin + ADVdin; % Bottom remineralization (ADVdin) added (Le Gland, 11/12/2019)
pondot = dPONdt + DIFFpon + ADVpon; %Only PON has vertical sinking.
%...................................................................................
boxdot = dBOXdt; 
%...................................................................................
%===================================================================================
if strcmp(keyPhysics,'not')
    %...............................................................................
    yavedot  = dYAVEdt;
    yyvardot = dYYVARdt;
    %...............................................................................
% $$$     xavedot = dXAVEdt + DIFFxave;
% $$$     xvardot = dXVARdt + DIFFxvar;
% $$$     xstddot = dXSTDdt + DIFFxstd;
    %...............................................................................
% $$$     xavedot = dXAVEdt + DIFFxaveout(:,jcounter);
% $$$     xvardot = dXVARdt + DIFFxvarout(:,jcounter);
% $$$     xstddot = dXSTDdt + DIFFxstdout(:,jcounter);
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    yave_stardot  = dYAVE_STARdt  + DIFFyave_star;
    yyvar_stardot = dYYVAR_STARdt + DIFFyyvar_star;
    %...............................................................................
end
%...................................................................................
%===================================================================================
%...................................................................................

%%%%%%%%%%
%STOCKAGE:
%%%%%%%%%%
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
UYday(:,jday) = uy;
GYday(:,jday) = gy;
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
%...................................................................................

%%%%%%%%
%OUTPUT:
%%%%%%%%
%===================================================================================
%...................................................................................
if strcmp(keyPhysics,'not')
    %...............................................................................
    Vdot = [phydot;zoodot;dindot;pondot;boxdot;yavedot;yyvardot];
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    Vdot = [phydot;zoodot;dindot;pondot;boxdot;yave_stardot;yyvar_stardot];
    %...............................................................................
end

return
%...................................................................................
%===================================================================================
%***********************************************************************************
