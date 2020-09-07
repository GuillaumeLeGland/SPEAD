function [Vdot] = jamstecrest_gaussecomodel1D_ode45eqs(iTime,V0)
global galfa gbeta
% global gzmax kgz mz betaz betap mpower 
global gzmax kgz mz betaz betap_max betap_xmax betap_xrange mpower 
global mp Isat InhFac numutx numuty % Replace numut by 2 different mutation rates (Le Gland, 10/09/2019) 
global alp0 mup0 % knp0 
global amup aknp bmup % aalp 
global amuz % Zooplankton dependent grazing (Le Gland, 26/08/2020)
%global Q10 
global Q10a Q10h % Distinct partition coefficients for auto and heterotrophic processes (Le Gland, 31/10/2019)
global ntot0 
global temp0 temp % sst % Use temperature at all depths (Le Gland, 15/10/2019) 
global jcounter 
%global iTimeode 
global keyPhysics keyAssimConstant % keyTraitAxis keySinking % keyAssimConstant added by Le Gland, 03/10/2019
global Sdin Drate % Mdin 
global epsPhy omePhy epsZoo omeZoo md 
global deltat % t0 ndays nyear tmax tspan 
global zdepths ndepths deltaz
global jday % jjday
%global Iave Ivar Istd %continuous model
global Ixave Iyave Ixxvar Iyyvar Ixycov % continuous model with 2 traits (Le Gland, 02/07/2019)
global Iphy Izoo Idin Ipon Ibox %continuous model
global KZ % KZI
global parz0
global kw wsink %kp
% global keyNutrientSupply 
%...................................................................................
% global Xavedotday Xvardotday Xstddotday %OUTPUTS 
% global UXday GXday 
% global d1UXdxday d1GXdxday 
% global d2UXdxday d2GXdxday 
global todedotday 
%...................................................................................
% global Xavedotout Xvardotout Xstddotout %OUTPUTS 
% global UXout GXout 
% global d1UXdxout d1GXdxout 
% global d2UXdxout d2GXdxout 
global todedotout 
% Outputs in the case with 2 traits (Le Gland, 16/07/2019)
% global Xavedotday Yavedotday XXvardotday YYvardotday XYcovdotday %OUTPUTS 
% global UXYday GXYday 
% global d1UXYdxday d1UXYdyday d1GXYdxday d1GXYdyday 
% global d2UXYdxdxday d2UXYdydyday d2UXYdxdyday d2GXYdxdxday d2GXYdydyday d2GXYdxdyday
% % Outputs in the case with 2 traits (Le Gland, 16/07/2019)
% global Xavedotout Yavedotout XXvardotout YYvardotout XYcovdotout %OUTPUTS 
global UXYout GXYout 
% global d1UXYdxout d1UXYdyout d1GXYdxout d1GXYdyout 
% global d2UXYdxdxout d2UXYdydyout d2UXYdxdyout d2GXYdxdxout d2GXYdydyout d2GXYdxdyout 
%...................................................................................
global FPHYToutcont EPHYToutcont MPHYToutcont GPHYToutcont %OUTPUTS 
global FZOOoutcont  EZOOoutcont  MZOOoutcont 
global FDINoutcont  FPONoutcont
%...................................................................................
%global DIFFxaveout DIFFxvarout DIFFxstdout 
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
    % $$$ xave = V0(Iave); %mean size
    % $$$ xvar = V0(Ivar); %variance
    % $$$ xstd = V0(Istd); %deviation
    % Case with 2 traits (Le Gland, 02/07/2019)
    xave  = V0(Ixave) ;
    yave  = V0(Iyave) ;
    xxvar = V0(Ixxvar);
    yyvar = V0(Iyyvar);
    xycov = V0(Ixycov);
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    % $$$ xave_star = V0(Iave); %mean size times phy.
    % $$$ xvar_star = V0(Ivar); %mean size squared plus variance times phy.
    % $$$ xstd_star = V0(Istd); %mean size plus standardeviation times phy.
    % Case with 2 traits (Le Gland, 02/07/2019)
    xave_star  = V0(Ixave) ;
    yave_star  = V0(Iyave) ;
    xxvar_star = V0(Ixxvar);
    yyvar_star = V0(Iyyvar);
    xycov_star = V0(Ixycov);
    %
    %...............................................................................
    %WRONG!!!
% $$$     xave = (xave_star./phy); %mean size
% $$$     xvar = (xvar_star./phy); %variance
% $$$     xstd = (xstd_star./phy); %deviation
    %...............................................................................
    %OKAY: 
    % xave = (xave_star./phy); %mean size
    % xvar = (xvar_star./phy) - xave.^2; %variance
    % xstd = (xstd_star./phy) - xave; %deviation
    % Case with 2 traits (Le Gland, 02/07/2019)
    xave  = (xave_star ./phy);
    yave  = (yave_star ./phy);
    xxvar = (xxvar_star./phy) - xave.^2   ;
    yyvar = (yyvar_star./phy) - yave.^2   ;
    xycov = (xycov_star./phy) - xave.*yave;
    %...............................................................................
end
% Protection against negative variances and out-of-range correlations (Le Gland, 21/11/2019)
varmin = 10^(-6); % To avoid eps*eps=0;
xxvar = max(varmin, xxvar); 
yyvar = max(varmin, yyvar);
%xycov = min(0.99*sqrt(xxvar.*yyvar), abs(xycov)).*sign(xycov);
% Avoid correlation smaller than 1 and add correction term (Le Gland, 26/02/2020)
dXYCOVdt_control = zeros(ndepths,1);
for i=1:ndepths
    %if sqrt(xxvar(i)*yyvar(i))-eps < abs(xycov(i))
    if sqrt(xxvar(i)*yyvar(i))*0.999 < abs(xycov(i))    % Stricter conditions to avoid bugs in the results (Le Gland, 20/07/2020)
        %dXYCOVdt_control(i) = (sqrt(xxvar(i)*yyvar(i)) - eps - abs(xycov(i)))*sign(xycov(i))/deltat;
        dXYCOVdt_control(i) = (sqrt(xxvar(i)*yyvar(i))*0.998 - abs(xycov(i)))*sign(xycov(i))/deltat;
        %xycov(i) = (sqrt(xxvar(i)*yyvar(i))-eps).*sign(xycov(i));
        xycov(i) = 0.999*sqrt(xxvar(i)*yyvar(i)).*sign(xycov(i));
    end
end


%...................................................................................
%===================================================================================
%................................................................................... 
VNPZD0 = [phy;zoo;din;pon]; 
%...................................................................................
%VSTAT0 = [xave;xvar;xstd];
%VSTAT0 = [xave,yave,xxvar,yyvar,xycov]; % Case with 2 traits (Le Gland, 02/07/2019)
%................................................................................... 
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DAY OF SIMULATION ANT TIME COUNTER:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
% jjday is identical to jday in all cases ! (Le Gland, 08/07/2019)
%[jday,jjday,newday] = jamstecrest_daycounter(jday,iTime);
[jday,newday] = jamstecrest_daycounter(jday,iTime); % I simplify the function (Le Gland, 13/09/2019)
%...................................................................................
%jcounter1 = jcounter; %Previou-s jcounter.
%...................................................................................
jcounter = floor(iTime/deltat); %For ode4.
% $$$ jcounter = floor(iTime/deltat) + 1; %For ode1.
%...................................................................................
%jcounter2 = jcounter; %Current jcounter.
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
xj = xave; 
xmj = xave; 
%...................................................................................
%%sigmaxj = xstd; %THIS I DONT REALLY KNOW WHY IS NOT WORKING.
%%sigmaxj = sqrt(xvar); %THIS SHOULD BE CORRECT ONE FOR SURE.
sigmaxj = sqrt(xxvar); % Le Gland, 05/07/2019

% Second trait (Le Gland, 02/07/2019)
yj = yave;
ymj = yave;
sigmayj = sqrt(yyvar);

%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK FOR NEGATIVE VARIANCE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%Ineg = find(xvar < 0);
%...................................................................................
%xvar(Ineg) = sqrt(eps); %To avoid negative variances -- Force them to be close to zero value.
%...................................................................................
%Jneg = find(xvar < 0);
%...................................................................................
% Case with 2 traits (Le Gland, 05/07/2019)
%Inegx = find(xxvar < 0);
%Inegy = find(yyvar < 0);
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
% If variance is zero, singular values are obtained (Le Gland, 12/06/2019)
% Since these derivatives are NOT necessary, they can be removed to avoid failure
% IDEA: instead of having a whole new mode for them, why not imposing a
% very small variance instead, like 0.001 or 10-9 ... (Le Gland, 18/09/2019)
% if sigmaxj > 0
    %fxj = (1.0 ./ (sigmaxj * sqrt(2*pi))) .* exp( -(xj - xmj).^2 ./ (2*sigmaxj.^2) ); %Okay.
    
    % Case for 2 traits, but it is just a ndepths*1 vector repeating 
    % the maximum value of the normal distribution (Le Gland, 02/07/2019)
    
% Simplify the case without KTW and avoid collapse when sigma is a singular 
% matrix (Le Gland, 26/02/2020)
tempj = temp(:,jday); % (Le Gland, 04/11/2019)
if sum(galfa ~= 1) ~= 0
        
    fxyj         = zeros(ndepths,1)  ;
    sigma        = zeros(ndepths,2,2);
    invsigma     = zeros(ndepths,2,2);
    detsigma     = zeros(ndepths,1)  ;
    detsigmasqrt = zeros(ndepths,1)  ;
    vecj         = zeros(ndepths,2)  ;
    for i=1:ndepths
       % Perhaps size and Topt should be expressed with two different trait units
       % But some matrices would have terms with different dimensions
       % Topt could be expressed as a temperature, since it can be compared
       % with environment temperatures, but plankton size must not be
       % expressed as a distance.
       sigma(i,:,:)        = [xxvar(i), xycov(i); xycov(i), yyvar(i)];  % Variance matrix [trait^2]
       % This is only necessary in KTW mode (Le Gland, 26/11/2019)
       % invsigma(i,:,:)     = inv(squeeze(sigma(i,:,:)));                % Inverse of variance matrix [trait^(-2)]
       detsigma(i)         = (xxvar(i) .* yyvar(i) - xycov(i).^2);      % Determinant of variance matrix [trait^4]
       detsigmasqrt(i)     = sqrt(detsigma(i));                         % Square root of variance determinant [trait^2]
       vecj(i,:)           = [xj(i)-xmj(i),yj(i)-ymj(i)];               % Trait vector compared to mean trait [trait]
       %vecjnow = squeeze(vecj(i,:));
    
       % fxyj(i) = (1.0 ./ (detsigmasqrt(i)*2*pi)) .* exp( -(1/2)*vecjnow*invsigma(i)*vecjnow' ); % 2D probability density function [trait^(-2)]
       % use \ instead of inverse matrix, to increase precision and decrease cost (Le Gland, 21/11/2019)
       % fxyj(i) = (1.0 ./ (detsigmasqrt(i)*2*pi)) .* exp( -(1/2)*vecjnow*((squeeze(sigma(i,:,:)\vecjnow') );
       % Anyway, this should be zero (!) (Le Gland, 21/11/2019)
       fxyj(i) = (1.0 ./ (detsigmasqrt(i)*2*pi));
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PHYTOPLANKTON BIOMASS GAUSSIAN DISTRIBUTION:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ptot = phy;
    %Pxj = Ptot .* fxj;
    %Pxjalfa = Pxj.^galfa;
    % Case with 2 traits (Le Gland, 03/07/2019)
    % Perhaps there should be an option for Kill The Winner on size only
    Pxyj = Ptot .* fxyj;
    Pxyjalfa = Pxyj.^galfa;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ANALYTICAL INTEGRATION OF PHYTOPLANKTON GASSIAN DISTRUBUTION:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===================================================================================
    %...................................................................................
    % intPxdx = Ptot;
    % intPxydxdy = Ptot; % 2 traits (Le Gland, 03/07/2019)
    %...................................................................................

    % intPxalfadx = (Ptot.^galfa ./ sqrt(galfa)) .* (sigmaxj*sqrt(2*pi)).^(1-galfa);
    intPxyalfadxdy = (Ptot.^galfa ./ galfa) .* (detsigmasqrt*2*pi).^(1-galfa); % 2 traits (Le Gland, 03/07/2019) 
    % intPxyalfadxdy = (Ptot.^galfa ./ sqrt(galfa)) .* (sigmaxj*sqrt(2*pi)).^(1-galfa); % KTW on size only (Le Gland, 27/09/2019)
    % intPxyalfadxdy = (Ptot.^galfa ./ sqrt(galfa)) .* (sigmayj*sqrt(2*pi)).^(1-galfa); % KTW on temperature only
    
    %...................................................................................
    %===================================================================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ANALYTICAL DERIVATIVES OF PHYTOPLANKTON GAUSSIAN DISTRIBUTION:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===================================================================================
    %...................................................................................
    % dPxdx = -(xj - xmj) .* (Pxj ./ sigmaxj.^2);
    %...................................................................................
    % d2Pxdx = Pxj .* ( ((xj - xmj).^2 ./ sigmaxj.^4) - (1.0 ./ sigmaxj.^2) );
    
    % 2 traits (Le Gland, 03/07/2019)
%     dPxydx    = zeros(ndepths,1);
%     dPxydy    = zeros(ndepths,1);
%     d2Pxydxdx = zeros(ndepths,1);
%     d2Pxydydy = zeros(ndepths,1);
%     d2Pxydxdy = zeros(ndepths,1);
%     for i=1:ndepths
%        dPxydx(i) = - Pxyj(i) * vecj(i,:) * squeeze(invsigma(i,1,:)) ; % Dimensions are dubious in this matrix product
%        dPxydy(i) = - Pxyj(i) * vecj(i,:) * squeeze(invsigma(i,2,:)) ; % to be improved (Le Gland, 10/07/2019)
%        %................................................................................
%        d2Pxydxdx(i) = Pxyj(i) * ( ( vecj(i,:) * squeeze(invsigma(i,1,:)) )^2 - invsigma(i,1,1) ) ;
%        d2Pxydydy(i) = Pxyj(i) * ( ( vecj(i,:) * squeeze(invsigma(i,2,:)) )^2 - invsigma(i,2,2) ) ;
%        d2Pxydxdy(i) = Pxyj(i) * ( ( vecj(i,:) * squeeze(invsigma(i,1,:)) )*(vecj(i,:) * squeeze(invsigma(i,2,:)) ) - invsigma(i,1,2) ) ;
%     end
    
    %...................................................................................
    %===================================================================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ANALYTICAL DERIVATIVES OF PHYTOPLANKTON-ALFA GAUSSIAN DISTRIBUTION:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===================================================================================
    %...................................................................................
    % dPxalfadx = galfa .* (Pxj.^galfa) .* (- (xj - xmj) ./ sigmaxj.^2); %Okay!!!
    %...................................................................................
    % d2Pxalfadx = galfa.*Pxjalfa .* (galfa.*((xj - xmj).^2 ./ sigmaxj.^4) - (1.0./sigmaxj.^2)); 
    %...................................................................................
    % 2 traits (Le Gland, 10/07/2019)
    % My derivation is different from the previous one
    % Still needs to be checked (identical to the previous formula in 1D case ?)
%     dPxyalfadx = - galfa .* (Pxyj.^(galfa-1)) .* dPxydx;
%     dPxyalfady = - galfa .* (Pxyj.^(galfa-1)) .* dPxydy;
    %...................................................................................
    % d2Pxyalfadxdx = galfa .* (galfa-1) .* (Pxyj.^(galfa-2)) .* dPxydx.^2      + galfa .* (Pxyj.^(galfa-1)) .* d2Pxydxdx;
    % d2Pxyalfadydy = galfa .* (galfa-1) .* (Pxyj.^(galfa-2)) .* dPxydy.^2      + galfa .* (Pxyj.^(galfa-1)) .* d2Pxydydy;
    % d2Pxyalfadxdy = galfa .* (galfa-1) .* (Pxyj.^(galfa-2)) .* dPxydx.*dPxydy + galfa .* (Pxyj.^(galfa-1)) .* d2Pxydxdy;
    %...................................................................................
    % The old derivation was the right one. It has to be reinstated in order to have the continuous model correspond to the
    % discrete model. Only the second order derivatives are different (Le Gland, 25/09/2019).
    % d2Pxyalfadxdx = galfa .* (Pxyj.^(galfa-1)) .* d2Pxydxdx;
    % d2Pxyalfadydy = galfa .* (Pxyj.^(galfa-1)) .* d2Pxydydy;
    % d2Pxyalfadxdy = galfa .* (Pxyj.^(galfa-1)) .* d2Pxydxdy;
%     for i=1:ndepths
%        d2Pxyalfadxdx(i) = galfa(i) * Pxyjalfa(i) * ( ( galfa(i) * vecj(i,:) * squeeze(invsigma(i,1,:)) )^2 - invsigma(i,1,1) ) ;
%        d2Pxyalfadydy(i) = galfa(i) * Pxyjalfa(i) * ( ( galfa(i) * vecj(i,:) * squeeze(invsigma(i,2,:)) )^2 - invsigma(i,2,2) ) ;
%        d2Pxyalfadxdy(i) = galfa(i) * Pxyjalfa(i) * ( ( galfa(i) * vecj(i,:) * squeeze(invsigma(i,1,:)) )*(vecj(i,:) * squeeze(invsigma(i,2,:)) ) - invsigma(i,1,2) ) ;
%     end
    %===================================================================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %GRAZING FUNCTIONAL RESPONSE KTW:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===================================================================================
    %...................................................................................
    % Vmax = (gzmax*zoo);
    % Add temperature effect (Le Gland, 21/10/2019)
    % Q10h for zooplankton is larger than Q10a for phytoplankton
    % Vmax = (gzmax*zoo).*(Q10h.^((tempj - temp0)/10));
    % Size dependent grazing (26/08/2020)
    Vmax = (gzmax*zoo) .* (Q10h.^((tempj - temp0)/10)) .* exp(amuz.*(xmj/10))
    % Qswitchj = (Pxjalfa ./ intPxalfadx); 
    % Qfeeding = (intPxdx.^gbeta ./ (intPxdx.^gbeta + kgz^gbeta)); %Okay.
    % 2 traits (Le Gland, 11/07/2019)
    Qswitchj = (Pxyjalfa ./ intPxyalfadxdy);
    %Qfeeding = (intPxydxdy.^gbeta ./ (intPxydxdy.^gbeta + kgz^gbeta));
    Qfeeding = (Ptot.^gbeta ./ (Ptot.^gbeta + kgz^gbeta));
    %...................................................................................
    % Gxj = Qswitchj .* Qfeeding .* Vmax;
    % 2 traits (Le Gland, 11/07/2019)
    Gxyj = Qswitchj .* Qfeeding .* Vmax;
    %===================================================================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ANALYTICAL DERIVATIVES OF GRAZING KTW FUNCTION:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===================================================================================
    %...................................................................................
    % A = (sigmaxj*sqrt(2*pi)).^(1-galfa) .* (1./sqrt(galfa));
    % 2 traits (Le Gland, 11/07/2019)
    % A = (sigmaxj.*sigmayj.*2.*pi).^(1-galfa) .* (1./sqrt(galfa));
    % A = (detsigmasqrt*2*pi).^(1-galfa) .* (1./galfa); % Correction by Le Gland (25/09/2019)
    % A = (sigmaxj*sqrt(2*pi)).^(1-galfa) .* (1./sqrt(galfa)); % KTW on size only (Le Gland, 27/09/2019)
    % A = (sigmayj*sqrt(2*pi)).^(1-galfa) .* (1./sqrt(galfa));
    % size only
    %...................................................................................
    %%G = (Ptot^gbeta / (Ptot^gbeta + kgz^gbeta)) * Vmax; 
    % G = (Qfeeding .* Vmax); 
    %...................................................................................
    % B = (G./A); 
    %...................................................................................
    % d1Gxdx = galfa .* fxj.^galfa .* (- (xj - xmj) ./ sigmaxj.^2) .* B;
    %...................................................................................
    % d2Gxdx = galfa .* fxj.^galfa .* (galfa .* ((xj - xmj).^2 ./ sigmaxj.^4) - (1./sigmaxj.^2)) .* B;
    % 2 traits (Le Gland, 11/07/2019)
%     d1Gxydx   = zeros(ndepths,1);
%     d1Gxydy   = zeros(ndepths,1);
%     d2Gxydxdx = zeros(ndepths,1);
%     d2Gxydydy = zeros(ndepths,1);
%     d2Gxydxdy = zeros(ndepths,1);
%     for i = 1:ndepths
%        d1Gxydx(i)   = - galfa(i) .* fxyj(i)^galfa(i) * (vecj(i,:) * squeeze(invsigma(i,1,:)) ) * B(i);
%        d1Gxydy(i)   = - galfa(i) .* fxyj(i)^galfa(i) * (vecj(i,:) * squeeze(invsigma(i,2,:)) ) * B(i);
%        d2Gxydxdx(i) = galfa(i) .* fxyj(i)^galfa(i) * ( galfa(i) * ( vecj(i,:) * squeeze(invsigma(i,1,:)) )^2 - invsigma(i,1,1) ) * B(i);
%        d2Gxydydy(i) = galfa(i) .* fxyj(i)^galfa(i) * ( galfa(i) * ( vecj(i,:) * squeeze(invsigma(i,2,:)) )^2 - invsigma(i,2,2) ) * B(i);
%        d2Gxydxdy(i) = galfa(i) .* fxyj(i)^galfa(i) * ( galfa(i) * ( vecj(i,:) * squeeze(invsigma(i,1,:)) )*(vecj(i,:) * squeeze(invsigma(i,2,:)) ) - invsigma(i,1,2) ) * B(i);
%     end
    %...................................................................................
    %===================================================================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %BIOMASS SPECIFIC GRAZING FUNCTIONAL RESPONSE KTW:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %...................................................................................
    % gxj = Gxj./Pxj; %[d-1] 
    % 2 traits (Le Gland, 11/07/2019)
    gxyj = Gxyj./Pxyj;
    %...................................................................................

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ANALYTICAL DERIVATIVES OF BIOMASS SPECIFIC GRAZING FUNCTIONAL RESPONSE KTW:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===================================================================================
    %...................................................................................
    % d1gxdx = gxj .* (galfa - 1) .* (-(xj - xmj) ./ sigmaxj.^2 ); 
    %...................................................................................
    % d2gxdx = gxj .* (galfa - 1) .* ( (galfa - 1) .* (-(xj - xmj) ./ sigmaxj.^2).^2 - (1./sigmaxj.^2) );
    %...................................................................................
    % 2 traits (Le Gland, 11/07/2019)
        
    d1gxydx   = zeros(ndepths,1);
    d1gxydy   = zeros(ndepths,1);
    d2gxydxdx = zeros(ndepths,1);
    d2gxydydy = zeros(ndepths,1);
    d2gxydxdy = zeros(ndepths,1);
    for i = 1:ndepths
       % d1gxydx(i)   = - gxyj(i) * (galfa(i) - 1) * (vecj(i,:) * squeeze(invsigma(i,1,:)) );
       d1gxydy(i)   = - gxyj(i) * (galfa(i) - 1) * (vecj(i,:) * squeeze(invsigma(i,2,:)) );
       % d2gxydxdx(i) = gxyj(i) * (galfa(i) - 1) * ( (galfa(i) - 1) * ( vecj(i,:) * squeeze(invsigma(i,1,:)) )^2 - invsigma(i,1,1) );
       d2gxydydy(i) = gxyj(i) * (galfa(i) - 1) * ( (galfa(i) - 1) * ( vecj(i,:) * squeeze(invsigma(i,2,:)) )^2 - invsigma(i,2,2) );
       d2gxydxdy(i) = gxyj(i) * (galfa(i) - 1) * ( (galfa(i) - 1) * ( vecj(i,:) * squeeze(invsigma(i,1,:)) )*(vecj(i,:) * squeeze(invsigma(i,2,:)) ) - invsigma(i,1,2) );
       % KTW on size only (Le Gland, 27/09/2019)
       % d2gxydxdx(i) = gxyj(i) * (galfa(i) - 1) * ( (galfa(i) - 1) * ( vecj(i,:) * squeeze(invsigma(i,1,:)) )^2 - invsigma(i,1,1) );
       % d2gxydydy(i) = 0;
       % d2gxydxdx(i) = 0;
       % d2gxydydy(i) = gxyj(i) * (galfa(i) - 1) * ( (galfa(i) - 1) * ( vecj(i,:) * squeeze(invsigma(i,2,:)) )^2 - invsigma(i,2,2) );
       % d2gxydxdy(i) = 0;
       % Added dependence on prey size (Le Gland, 26/08/2020)
       % d1gxdx = gxj .* (galfa - 1) .* (-(xj - xmj) ./ sigmaxj.^2 ) + (amuz .* gxj);
       % d2gxdx = gxj .* (galfa - 1) .* ( (galfa - 1) .* (-(xj - xmj) ./ sigmaxj.^2).^2 - (1./sigmaxj.^2) ) + ...
       %         2 .* amuz .* gxj .* (galfa - 1) .* (-(xj - xmj) ./ sigmaxj.^2 ) + (amuz^2 .* gxj);
       d1gxydx(i)   = - gxyj(i) .* (galfa(i) - 1) .* (vecj(i,:) * squeeze(invsigma(i,1,:)) ) + (amuz .* gxyj(i));
       d2gxydxdx(i) = gxyj(i) .* (galfa(i) - 1) .* ( (galfa(i) - 1) .* (vecj(i,:) * invsigma(i,1,:)).^2 - invsigma(i,1,1) ) + ...
                      2 .* amuz .* gxyj(i) .* (galfa - 1) .* (- vecj(i,:) * squeeze(invsigma(i,1,:)) ) + (amuz^2 .* gxyj(i));
    end
elseif sum(galfa ~= 1) == 0
    Vmax = (gzmax*zoo) .* (Q10h.^((tempj - temp0)/10)) .* exp(amuz.*(xmj/10));
    Qfeeding = (phy.^gbeta ./ (phy.^gbeta + kgz^gbeta));
    d1gxydx   = 0;
    d1gxydy   = 0;
    d2gxydxdx = 0;
    d2gxydydy = 0;
    d2gxydxdy = 0;
    % Added dependence on prey size (Le Gland, 14/10/2019)
    % d1gxdx = gxj .* (galfa - 1) .* (-(xj - xmj) ./ sigmaxj.^2 ) + (amuz .* gxj);
    % d2gxdx = gxj .* (galfa - 1) .* ( (galfa - 1) .* (-(xj - xmj) ./ sigmaxj.^2).^2 - (1./sigmaxj.^2) ) + ...
    %          2 .* amuz .* gxj .* (galfa - 1) .* (-(xj - xmj) ./ sigmaxj.^2 ) + (amuz^2 .* gxj);
    fxyj         = zeros(ndepths,1)  ;
    sigma        = zeros(ndepths,2,2);
    invsigma     = zeros(ndepths,2,2);
    detsigma     = zeros(ndepths,1)  ;
    detsigmasqrt = zeros(ndepths,1)  ;
    vecj         = zeros(ndepths,2)  ;
    for i=1:ndepths
       % Perhaps size and Topt should be expressed with two different trait units
       % But some matrices would have terms with different dimensions
       % Topt could be expressed as a temperature, since it can be compared
       % with environment temperatures, but plankton size must not be
       % expressed as a distance.
       sigma(i,:,:)        = [xxvar(i), xycov(i); xycov(i), yyvar(i)];  % Variance matrix [trait^2]
       % This is only necessary in KTW mode (Le Gland, 26/11/2019)
       % invsigma(i,:,:)     = inv(squeeze(sigma(i,:,:)));                % Inverse of variance matrix [trait^(-2)]
       detsigma(i)         = (xxvar(i) .* yyvar(i) - xycov(i).^2);      % Determinant of variance matrix [trait^4]
       detsigmasqrt(i)     = sqrt(detsigma(i));                         % Square root of variance determinant [trait^2]
       vecj(i,:)           = [xj(i)-xmj(i),yj(i)-ymj(i)];               % Trait vector compared to mean trait [trait]
       %vecjnow = squeeze(vecj(i,:));
       % fxyj(i) = (1.0 ./ (detsigmasqrt(i)*2*pi)) .* exp( -(1/2)*vecjnow*invsigma(i)*vecjnow' ); % 2D probability density function [trait^(-2)]
       % use \ instead of inverse matrix, to increase precision and decrease cost (Le Gland, 21/11/2019)
       % fxyj(i) = (1.0 ./ (detsigmasqrt(i)*2*pi)) .* exp( -(1/2)*vecjnow*((squeeze(sigma(i,:,:)\vecjnow') );
       % Anyway, this should be zero (!) (Le Gland, 21/11/2019)
       fxyj(i) = (1.0 ./ (detsigmasqrt(i)*2*pi));
    end
    %detsigma = xxvar.*yyvar-xycov.^2;
    %fxyj = (1.0 ./ (sqrt(detsigma)*2*pi));
    Qswitchj = fxyj;
    Gxyj = Qswitchj .* Qfeeding .* Vmax;
    Pxyj = phy.*fxyj;
    gxyj = Gxyj./Pxyj;
    d1gxydx   = amuz .* gxyj;
    d2gxydxdx = amuz^2 .* gxyj;
end
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
DIFFxave_star  = zeros(ndepths,1);
DIFFxxvar_star = zeros(ndepths,1);
%DIFFxstd_star = zeros(ndepths,1);
% Case with 2 traits (Le Gland, 03/07/2019)
DIFFyave_star  = zeros(ndepths,1);
DIFFyyvar_star = zeros(ndepths,1);
DIFFxycov_star = zeros(ndepths,1);
%...................................................................................
if ndepths > 1
    % [DIFFxave_star]  = jamstecrest_TurbulentDiffusion(xave_star,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
    % [DIFFxxvar_star] = jamstecrest_TurbulentDiffusion(xxvar_star,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
    %%[DIFFxstd_star] = jamstecrest_TurbulentDiffusion(xstd_star,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
    % Case with 2 traits (Le Gland, 03/07/2019)
    % [DIFFyave_star]  = jamstecrest_TurbulentDiffusion(yave_star,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
    % [DIFFyyvar_star] = jamstecrest_TurbulentDiffusion(yyvar_star,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
    % [DIFFxycov_star] = jamstecrest_TurbulentDiffusion(xycov_star,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
    %disp(DIFFxxvar_star)
    % kzI is no longer required (Le Gland, 10/09/2019)
    [DIFFxave_star]  = jamstecrest_TurbulentDiffusion(xave_star,deltat,deltaz,kz,ndepths,'Reflectante');
    [DIFFxxvar_star] = jamstecrest_TurbulentDiffusion(xxvar_star,deltat,deltaz,kz,ndepths,'Reflectante');
    [DIFFyave_star]  = jamstecrest_TurbulentDiffusion(yave_star,deltat,deltaz,kz,ndepths,'Reflectante');
    [DIFFyyvar_star] = jamstecrest_TurbulentDiffusion(yyvar_star,deltat,deltaz,kz,ndepths,'Reflectante');
    [DIFFxycov_star] = jamstecrest_TurbulentDiffusion(xycov_star,deltat,deltaz,kz,ndepths,'Reflectante');
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
    % [DIFFphy] = jamstecrest_TurbulentDiffusion(phy,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
    % [DIFFzoo] = jamstecrest_TurbulentDiffusion(zoo,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
    % [DIFFdin] = jamstecrest_TurbulentDiffusion(din,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
    % [DIFFpon] = jamstecrest_TurbulentDiffusion(pon,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
    % kzI is no longer required (Le Gland, 10/09/2019)
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
    %....................................................................
    % ADVxave  = zeros(ndepths,1);
    % ADVxxvar = zeros(ndepths,1);
    % ADVxstd = zeros(ndepths,1);
    % Case with 2 traits (Le Gland, 03/07/2019)
    % ADVyave  = zeros(ndepths,1);
    % ADVyyvar = zeros(ndepths,1);
    % ADVxycov = zeros(ndepths,1);
    %....................................................................
% elseif strcmp(keyPhysics,'yes')
    %....................................................................
    % ADVxave_star  = zeros(ndepths,1);
    % ADVxxvar_star = zeros(ndepths,1);
    % ADVxstd_star = zeros(ndepths,1);
    % Case with 2 traits (Le Gland, 03/07/2019)
    % ADVyave_star  = zeros(ndepths,1);
    % ADVyyvar_star = zeros(ndepths,1);
    % ADVxycov_star = zeros(ndepths,1);
    %....................................................................
% end
%........................................................................
%========================================================================
%........................................................................
% ADVphy = zeros(ndepths,1);
% ADVzoo = zeros(ndepths,1);
% ADVdin = zeros(ndepths,1);
% ADVpon = zeros(ndepths,1);
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
% Article way (24/08/2020)
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

% Case with both traits activated (Le Gland, 11/07/2019)

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
%knxj = exp(xj); % Case where trait is half-saturation (Le Gland, 26/06/2019)
%...............................................................................
%===============================================================================
%EXPONENTIAL MAXIMUM GROWTH RATE: 
%...............................................................................
% $$$     muxj = mup .* exp(amup*xj); %Phy maximum grazing rate as a function of cell size [d-1]
%...............................................................................
%===============================================================================
%UNIMODAL MAXIMUM GROWTH RATE AT REFERENCE TEMPERATURE (15 DEGREES):
%...............................................................................
muxj = mup .* exp(amup*xj + bmup*xj.^2);
%muxj = mup .* exp(amup*(xj-log(knp))/aknp); % Case where trait is half-saturation (Le Gland, 26/06/2019)
%...............................................................................
%===============================================================================
%...............................................................................
lxj = (knxj ./ (knxj + din)); %[n.d.]
qxj = (din  ./ (din + knxj)); %[n.d.]
%...............................................................................
% gammay = 4.0; %Tolerance range [Celsius]
% $$$ %%gammay = 1d6; %Control case without gaussian.
%............................................................................... 
%sstj = sst(jday); 
%sstj = temp(:,jday); % Use temperature at all depths (Le Gland, 15/10/2019)
%...............................................................................
%q10 = Q10.^((sstj - sst0)/10); 
%...............................................................................
% qyj = exp(-(xj - sstj).^2 / (2*gammay^2)) * q10; %SST limitation [%]
% A mistake in qyj produced spurious results (too low growth ant Topt tending to logsize values)
% (Le Gland, 19/07/2019)
% qyj = exp(-(yj - sstj).^2 / (2*gammay^2)) .* q10;


% New skewed response to temperature (Le Gland, 23/10/2019)
% q10 = Q10.^((yj - sst0)/10);
q10 = Q10a.^((yj - temp0)/10); % (Le Gland, 04/11/2019)

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
%d1qxdx = - qxj .* lxj; % Case where trait is half-saturation (Le Gland, 26/06/2019)
d1lxdx = - d1qxdx; 
%...............................................................................
d2qxdx = -aknp * (d1lxdx .* qxj + d1qxdx .* lxj); 
%d2qxdx = d1lxdx .* qxj + d1qxdx .* lxj; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%d2lxdx = - d2qxdx; 
%...............................................................................
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
cff = (amup + 2*bmup*xj); 
% cff = amup / aknp; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%...............................................................................
d1muxdx = muxj.*cff;
d2muxdx = muxj*2*bmup + muxj.*cff.^2;
%d3muxdx = (2*bmup+cff.^2).*d1muxdx + 4.*bmup*muxj.*cff;
%d4muxdx = d1muxdx*8*bmup.*cff + (2*bmup+cff.^2).*d2muxdx + 8*bmup.^2*muxj;
%d2muxdx = muxj.*cff.^2; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%d3muxdx = muxj.*cff.^3; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%d4muxdx = muxj.*cff.^4; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%...............................................................................
%===============================================================================
%FROM BINGZANG CHEN FOR NUTRIENT UPTAKE USING UNIMODAL MAXIMUM GROWTH RATE: 
%...............................................................................
d1uxydx = qyj .* ( (d1muxdx .* qxj) + (muxj .* d1qxdx) );
%...............................................................................
d2uxydxdx = qyj .* ( ((d2muxdx .*    qxj) + (d1muxdx .* d1qxdx)) + ...
          ((d1muxdx .* d1qxdx) + (muxj    .* d2qxdx)) );
%...............................................................................
%d3uxydx = qyj .* ( ((d3muxdx .*    qxj) + (d2muxdx .* d1qxdx)) + ...
%          ((d2muxdx .* d1qxdx) + (d1muxdx .* d2qxdx)) + ...
%          ((d2muxdx .* d1qxdx) + (d1muxdx .* d2qxdx)) + ...
%          ((d1muxdx .* d2qxdx) + (muxj    .* d3qxdx)) );
%...............................................................................
%d4uxydx = qyj .* ( ((d4muxdx .*    qxj) + (d3muxdx .* d1qxdx)) + ...
%          ((d3muxdx .* d1qxdx) + (d2muxdx .* d2qxdx)) + ...
%          ((d3muxdx .* d1qxdx) + (d2muxdx .* d2qxdx)) + ...
%          ((d2muxdx .* d2qxdx) + (d1muxdx .* d3qxdx)) + ...
%          ((d3muxdx .* d1qxdx) + (d2muxdx .* d2qxdx)) + ...
%          ((d2muxdx .* d2qxdx) + (d1muxdx .* d3qxdx)) + ...
%          ((d2muxdx .* d2qxdx) + (d1muxdx .* d3qxdx)) + ...
%          ((d1muxdx .* d3qxdx) + (muxj    .* d4qxdx)) );
%...............................................................................
%===============================================================================

% d1uxydy = uxyj .* (  -(yj - tempj)/gammay^2);
% d2uxydydy = uxyj .* ( (-(yj - tempj)/gammay^2).^2 - (1/gammay^2) );
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
    d1uxydy = uxyj .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 );
    d2uxydydy = d1uxydy .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 ) - uxyj .* 1./(yj + 5 - tempj).^2;    
    d2uxydxdy = qyj .* ( (d1muxdx .* qxj) + (muxj .* d1qxdx) ) .* ( - 0.2 + 1./(yj + 5 - tempj) + log(Q10a)/10 );
else
    d1uxydy = 0;
    d2uxydydy = 0;
    d2uxydxdy = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USING LAN SMITH ANALYTICAL DERIVATIVES: (JUST FOR CHECKING WITH HIM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================================
% $$$     %...............................................................................
% $$$     %%gxjLan = (Vmax * Qfeeding) * (1 / Ptot) * ( sqrt(galfa) * exp(-(1/2)*(galfa-1)*((xj-xmj)/sigmaxj)^2))
% $$$     gxjLan = (Vmax .* Qfeeding) .* (1 ./ Ptot) .* sqrt(galfa); 
% $$$     %...............................................................................
% $$$     d1gxdxLan = zeros(ndepths,1); 
% $$$     %...............................................................................
% $$$     d2gxdxLan = -(Vmax .* Qfeeding) .* (1 ./ Ptot) .* ( (galfa - 1) .* sqrt(galfa) ./ sigmaxj.^2 );
% $$$     %...............................................................................
% $$$     Qdin = qxj; %[n.d.]
% $$$     %...............................................................................
% $$$     uxjLan = muxj .* Qdin; 
% $$$     %...............................................................................
% $$$     d1uxdxLan = (amup*muxj.*Qdin) - aknp*(knxj./(muxj.*din)).*(muxj.*Qdin).^2;
% $$$     %...............................................................................
% $$$     d2uxdxLan = (amup - 2*aknp*((knxj.*muxj.*Qdin)./(muxj.*din))) .* d1uxdxLan ...
% $$$ 	- aknp*(aknp - amup)*(knxj./(muxj.*din)).*(muxj.*Qdin).^2;
% $$$     %...............................................................................
% $$$     d3uxdxLan = (amup - 2*aknp*((knxj.*muxj.*Qdin)./(muxj.*din))) .* d2uxdxLan ...
% $$$ 	- 2*aknp*(knxj./(muxj.*din)) .* d1uxdxLan.^2 ...
% $$$         - 4*aknp*(aknp - amup)    .* ((knxj.*muxj.*Qdin)./(muxj.*din)) .* d1uxdxLan ...
% $$$ 	-   aknp*(aknp - amup).^2 .* ( knxj./(muxj.*din)) .* uxjLan.^2;
% $$$     %...............................................................................
% $$$     d4uxdxLan = (amup - 2*aknp*(muxj.*Qdin.*knxj./(muxj.*din))) .* d3uxdxLan ...
% $$$         - 6*(aknp*(aknp-amup))*(muxj.*Qdin.*knxj./(muxj.*din)) .* d2uxdxLan ...
% $$$ 	- 6*(aknp*(knxj./(muxj.*din))) .* (d1uxdxLan .* d2uxdxLan) ...
% $$$ 	- 6*(aknp*(aknp-amup))*(knxj./(muxj.*din)) .* d1uxdxLan.^2 ...
% $$$ 	- 6*(aknp*(aknp-amup).^2)*(muxj.*Qdin.*knxj./(muxj.*din)) .* d1uxdxLan ...
% $$$ 	- (aknp*(aknp-amup).^3)*(knxj./(muxj.*din)) .* (muxj.*Qdin).^2;
% $$$     %...............................................................................
% $$$     %%gxj = gxjLan;
% $$$     %%d1gxdx = d1gxdxLan;
% $$$     %%d2gxdx = d2gxdxLan;
% $$$     %...............................................................................
% $$$     %%uxj = uxjLan;
% $$$     %%d1uxdx = d1uxdxLan;
% $$$     %%d2uxdx = d2uxdxLan;
% $$$     %...............................................................................
    %===============================================================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %USING MATLAB SYMBOLIC TOOLBOX DERIVATIVES: (JUST FOR DOUBLE CHECK)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===============================================================================
    %...............................................................................
% $$$     [d1uxdxSymbolic,d2uxdxSymbolic,d3uxdxSymbolic,d4uxdxSymbolic] = jamstecrest_gaussecomodel1D_symbolicderivatives;
% $$$     %...............................................................................
% $$$     d1uxdxSym = eval(d1uxdxSymbolic); %Convert from symbolic to numerical values at double precision.
% $$$     d2uxdxSym = eval(d2uxdxSymbolic);
% $$$     d3uxdxSym = eval(d3uxdxSymbolic);
% $$$     d4uxdxSym = eval(d4uxdxSymbolic);
    %..............................................................................
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %..............................................................................
% $$$     show_dudx1 = [d1uxdx,d1uxdxSym]
% $$$     show_dudx2 = [d2uxdx,d2uxdxSym]
% $$$     show_dudx3 = [d3uxdx,d3uxdxSym]
% $$$     show_dudx4 = [d4uxdx,d4uxdxSym]
% $$$     %..............................................................................
% $$$ % $$$     show_dudx1 = [d1uxdx,d1uxdxLan]
% $$$ % $$$     show_dudx2 = [d2uxdx,d2uxdxLan]
% $$$ % $$$     show_dudx3 = [d3uxdx,d3uxdxLan]
% $$$ % $$$     show_dudx4 = [d4uxdx,d4uxdxLan]
% $$$     %..............................................................................
% $$$ % $$$     show_dudx1 = [d1uxdx,d1uxdxLan,d1uxdxSym]
% $$$ % $$$     show_dudx2 = [d2uxdx,d2uxdxLan,d2uxdxSym]
% $$$ % $$$     show_dudx3 = [d3uxdx,d3uxdxLan,d3uxdxSym]
% $$$ % $$$     show_dudx4 = [d4uxdx,d4uxdxLan,d4uxdxSym]
% $$$     %..............................................................................
% $$$     disp('----')
% $$$     pause
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %===============================================================================

% if strcmp(keyTraitAxis,'ESD')
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %PHYTOPLANKTON GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %===============================================================================
%     %...............................................................................
%     knp = (mup./alp); 
%     %...............................................................................
%     %===============================================================================
%     %EXPONENTIAL UPTAKE AFFINITY: 
%     %...............................................................................
%     knxj = knp .* exp(aknp*xj); %Phy half-sat uptake as a function of cell size [mmolN*m-3]
%     %knxj = exp(xj); % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     %...............................................................................
%     %===============================================================================
%     %EXPONENTIAL MAXIMUM GROWTH RATE: 
%     %...............................................................................
% % $$$     muxj = mup .* exp(amup*xj); %Phy maximum grazing rate as a function of cell size [d-1]
%     %...............................................................................
%     %===============================================================================
%     %UNIMODAL MAXIMUM GROWTH RATE:
%     %...............................................................................
%     muxj = mup .* exp(amup*xj + bmup*xj.^2);
%     %muxj = mup .* exp(amup*(xj-log(knp))/aknp); % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     %...............................................................................
%     %===============================================================================
%     %...............................................................................
%     lxj = (knxj ./ (knxj + din)); %[n.d.]
%     qxj = (din  ./ (din + knxj)); %[n.d.]
%     %...............................................................................
%     uxj = muxj .* qxj; %[d-1] 
%     %...............................................................................
%     %===============================================================================
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %ANALYTICAL DERIVATIVES OF GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %===============================================================================
%     %MY DERIVATIONS FOR MICHAELIS MENTEN: 
%     %...............................................................................
%     d1qxdx = - qxj .* (aknp * lxj);
%     %d1qxdx = - qxj .* lxj; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     d1lxdx = - d1qxdx; 
%     %...............................................................................
%     d2qxdx = -aknp * (d1lxdx .* qxj + d1qxdx .* lxj); 
%     %d2qxdx = d1lxdx .* qxj + d1qxdx .* lxj; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     d2lxdx = - d2qxdx; 
%     %...............................................................................
%     d3qxdx = -aknp * (d2lxdx .* qxj + d2qxdx .* lxj + 2 * d1lxdx .* d1qxdx); 
%     %d3qxdx = d2lxdx .* qxj + d2qxdx .* lxj + 2 * d1lxdx .* d1qxdx; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     d3lxdx = - d3qxdx; 
%     %...............................................................................
%     d4qxdx = -aknp * ((d3lxdx .*   qxj  + d2lxdx .* d1qxdx) + ...
% 	              (d3qxdx .*   lxj  + d2qxdx .* d1lxdx) + ...
% 	              (d2lxdx .* d1qxdx + d1lxdx .* d2qxdx) * 2);
%     %d4qxdx = ((d3lxdx .*   qxj  + d2lxdx .* d1qxdx) + ...
% 	%              (d3qxdx .*   lxj  + d2qxdx .* d1lxdx) + ...
% 	%              (d2lxdx .* d1qxdx + d1lxdx .* d2qxdx) * 2); % Case where trait is half-saturation (Le Gland, 26/06/2019)
%               
%     %...............................................................................
%     %===============================================================================
%     %FROM BINGZANG CHEN FOR MICHAELIS MENTEN: 
%     %...............................................................................
% % $$$     d1qxdx_bis = -aknp*knxj.*qxj.^2;
% % $$$     d2qxdx_bis = -aknp^2*din.*knxj.*(1./(din+knxj).^2 - 2*knxj./(din+knxj).^3);
% % $$$     d3qxdx_bis =  aknp^3*din.*knxj.*(2*din.*knxj-(knxj-din).^2)./(knxj+din).^4;
% % $$$     d4qxdx_bis =  aknp^4*din.*knxj.*(11*knxj.*din.*(din-knxj)+knxj.^3-din.^3)./(din+knxj).^5; 
%     %...............................................................................
% % $$$     d1lxdx = - d1qxdx_bis; 
% % $$$     d2lxdx = - d2qxdx_bis; 
% % $$$     d3lxdx = - d3qxdx_bis; 
%     %...............................................................................
%     %===============================================================================
%     %MY DERIVATIONS FOR NUTRIENT UPTAKE USING EXPONENTIAL GROWTH RATE WITH SIZE: 
%     %...............................................................................
% % $$$     gff = (amup - aknp * lxj);
% % $$$     %...............................................................................
% % $$$     d1uxdx =   uxj  .*  gff; 
% % $$$     %...............................................................................
% % $$$     d2uxdx =   uxj  .* (gff.^2 - aknp^2 * lxj .* qxj); 
% % $$$     %...............................................................................
% % $$$     d3uxdx = d2uxdx .*  gff - 2*aknp * (d1uxdx .* d1lxdx) - aknp * (uxj .* d2lxdx); 
% % $$$     %...............................................................................
% % $$$     d4uxdx = d3uxdx .*  gff - 3*aknp * (d2uxdx .* d1lxdx + d1uxdx .* d2lxdx) - aknp * (uxj .* d3lxdx); 
%     %...............................................................................
%     %===============================================================================
%     %FROM BINGZANG CHEN FOR UNIMODAL MAXIMUM GROWTH RATE WITH SIZE: 
%     %...............................................................................
%     cff = (amup + 2*bmup*xj); 
%     % cff = amup / aknp; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     %...............................................................................
%     d1muxdx = muxj.*cff;
%     d2muxdx = muxj*2*bmup + muxj.*cff.^2;
%     d3muxdx = (2*bmup+cff.^2).*d1muxdx + 4.*bmup*muxj.*cff;
%     d4muxdx = d1muxdx*8*bmup.*cff + (2*bmup+cff.^2).*d2muxdx + 8*bmup.^2*muxj;
%     %d2muxdx = muxj.*cff.^2; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     %d3muxdx = muxj.*cff.^3; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     %d4muxdx = muxj.*cff.^4; % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     %...............................................................................
%     %===============================================================================
%     %FROM BINGZANG CHEN FOR NUTRIENT UPTAKE USING UNIMODAL MAXIMUM GROWTH RATE: 
%     %...............................................................................
%     d1uxdx = (d1muxdx .* qxj) + (muxj .* d1qxdx);
%     %...............................................................................
%     d2uxdx = ((d2muxdx .*    qxj) + (d1muxdx .* d1qxdx)) + ...
% 	     ((d1muxdx .* d1qxdx) + (muxj    .* d2qxdx));
%     %...............................................................................
%     d3uxdx = ((d3muxdx .*    qxj) + (d2muxdx .* d1qxdx)) + ...
% 	     ((d2muxdx .* d1qxdx) + (d1muxdx .* d2qxdx)) + ...
% 	     ((d2muxdx .* d1qxdx) + (d1muxdx .* d2qxdx)) + ...
% 	     ((d1muxdx .* d2qxdx) + (muxj    .* d3qxdx));
%     %...............................................................................
%     d4uxdx = ((d4muxdx .*    qxj) + (d3muxdx .* d1qxdx)) + ...
% 	     ((d3muxdx .* d1qxdx) + (d2muxdx .* d2qxdx)) + ...
% 	     ((d3muxdx .* d1qxdx) + (d2muxdx .* d2qxdx)) + ...
% 	     ((d2muxdx .* d2qxdx) + (d1muxdx .* d3qxdx)) + ...
% 	     ((d3muxdx .* d1qxdx) + (d2muxdx .* d2qxdx)) + ...
% 	     ((d2muxdx .* d2qxdx) + (d1muxdx .* d3qxdx)) + ...
% 	     ((d2muxdx .* d2qxdx) + (d1muxdx .* d3qxdx)) + ...
% 	     ((d1muxdx .* d3qxdx) + (muxj    .* d4qxdx));
%     %...............................................................................
%     %===============================================================================
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %USING LAN SMITH ANALYTICAL DERIVATIVES: (JUST FOR CHECKING WITH HIM)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %===============================================================================
% % $$$     %...............................................................................
% % $$$     %%gxjLan = (Vmax * Qfeeding) * (1 / Ptot) * ( sqrt(galfa) * exp(-(1/2)*(galfa-1)*((xj-xmj)/sigmaxj)^2))
% % $$$     gxjLan = (Vmax .* Qfeeding) .* (1 ./ Ptot) .* sqrt(galfa); 
% % $$$     %...............................................................................
% % $$$     d1gxdxLan = zeros(ndepths,1); 
% % $$$     %...............................................................................
% % $$$     d2gxdxLan = -(Vmax .* Qfeeding) .* (1 ./ Ptot) .* ( (galfa - 1) .* sqrt(galfa) ./ sigmaxj.^2 );
% % $$$     %...............................................................................
% % $$$     Qdin = qxj; %[n.d.]
% % $$$     %...............................................................................
% % $$$     uxjLan = muxj .* Qdin; 
% % $$$     %...............................................................................
% % $$$     d1uxdxLan = (amup*muxj.*Qdin) - aknp*(knxj./(muxj.*din)).*(muxj.*Qdin).^2;
% % $$$     %...............................................................................
% % $$$     d2uxdxLan = (amup - 2*aknp*((knxj.*muxj.*Qdin)./(muxj.*din))) .* d1uxdxLan ...
% % $$$ 	- aknp*(aknp - amup)*(knxj./(muxj.*din)).*(muxj.*Qdin).^2;
% % $$$     %...............................................................................
% % $$$     d3uxdxLan = (amup - 2*aknp*((knxj.*muxj.*Qdin)./(muxj.*din))) .* d2uxdxLan ...
% % $$$ 	- 2*aknp*(knxj./(muxj.*din)) .* d1uxdxLan.^2 ...
% % $$$         - 4*aknp*(aknp - amup)    .* ((knxj.*muxj.*Qdin)./(muxj.*din)) .* d1uxdxLan ...
% % $$$ 	-   aknp*(aknp - amup).^2 .* ( knxj./(muxj.*din)) .* uxjLan.^2;
% % $$$     %...............................................................................
% % $$$     d4uxdxLan = (amup - 2*aknp*(muxj.*Qdin.*knxj./(muxj.*din))) .* d3uxdxLan ...
% % $$$         - 6*(aknp*(aknp-amup))*(muxj.*Qdin.*knxj./(muxj.*din)) .* d2uxdxLan ...
% % $$$ 	- 6*(aknp*(knxj./(muxj.*din))) .* (d1uxdxLan .* d2uxdxLan) ...
% % $$$ 	- 6*(aknp*(aknp-amup))*(knxj./(muxj.*din)) .* d1uxdxLan.^2 ...
% % $$$ 	- 6*(aknp*(aknp-amup).^2)*(muxj.*Qdin.*knxj./(muxj.*din)) .* d1uxdxLan ...
% % $$$ 	- (aknp*(aknp-amup).^3)*(knxj./(muxj.*din)) .* (muxj.*Qdin).^2;
% % $$$     %...............................................................................
% % $$$     %%gxj = gxjLan;
% % $$$     %%d1gxdx = d1gxdxLan;
% % $$$     %%d2gxdx = d2gxdxLan;
% % $$$     %...............................................................................
% % $$$     %%uxj = uxjLan;
% % $$$     %%d1uxdx = d1uxdxLan;
% % $$$     %%d2uxdx = d2uxdxLan;
% % $$$     %...............................................................................
%     %===============================================================================
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %USING MATLAB SYMBOLIC TOOLBOX DERIVATIVES: (JUST FOR DOUBLE CHECK)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %===============================================================================
%     %...............................................................................
% % $$$     [d1uxdxSymbolic,d2uxdxSymbolic,d3uxdxSymbolic,d4uxdxSymbolic] = jamstecrest_gaussecomodel1D_symbolicderivatives;
% % $$$     %...............................................................................
% % $$$     d1uxdxSym = eval(d1uxdxSymbolic); %Convert from symbolic to numerical values at double precision.
% % $$$     d2uxdxSym = eval(d2uxdxSymbolic);
% % $$$     d3uxdxSym = eval(d3uxdxSymbolic);
% % $$$     d4uxdxSym = eval(d4uxdxSymbolic);
%     %..............................................................................
%     %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%     %..............................................................................
% % $$$     show_dudx1 = [d1uxdx,d1uxdxSym]
% % $$$     show_dudx2 = [d2uxdx,d2uxdxSym]
% % $$$     show_dudx3 = [d3uxdx,d3uxdxSym]
% % $$$     show_dudx4 = [d4uxdx,d4uxdxSym]
% % $$$     %..............................................................................
% % $$$ % $$$     show_dudx1 = [d1uxdx,d1uxdxLan]
% % $$$ % $$$     show_dudx2 = [d2uxdx,d2uxdxLan]
% % $$$ % $$$     show_dudx3 = [d3uxdx,d3uxdxLan]
% % $$$ % $$$     show_dudx4 = [d4uxdx,d4uxdxLan]
% % $$$     %..............................................................................
% % $$$ % $$$     show_dudx1 = [d1uxdx,d1uxdxLan,d1uxdxSym]
% % $$$ % $$$     show_dudx2 = [d2uxdx,d2uxdxLan,d2uxdxSym]
% % $$$ % $$$     show_dudx3 = [d3uxdx,d3uxdxLan,d3uxdxSym]
% % $$$ % $$$     show_dudx4 = [d4uxdx,d4uxdxLan,d4uxdxSym]
% % $$$     %..............................................................................
% % $$$     disp('----')
% % $$$     pause
%     %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%     %===============================================================================
% elseif strcmp(keyTraitAxis,'SST')
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %USING NEW TRAIT Y-AXIS FOR THE OPTIMAL TEMPERATURE FOR GROWTH:
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %===============================================================================
%     %...............................................................................
%     % yj = xj; 
%     % ymj = xmj; 
%     %...............................................................................
%     gammay = 4.0; %Tolerance range [Celsius]
% % $$$     %%gammay = 1d6; %Control case without gaussian.
%     %...............................................................................
% % $$$     %%sstj = sst(jcounter); 
%     sstj = sst(jday); 
%     %...............................................................................
%     q10 = Q10 ^((sstj - sst0)/10); 
%     %...............................................................................
%     qyj = exp(-(xj - sstj).^2 / (2*gammay^2)) * q10; %SST limitation [%]
%     %...............................................................................
%     knp = knp0*ones(ndepths,1);
%     %...............................................................................
%     uxj = mup .* qyj .* (din ./ (din + knp)); %Uptake rate [d-1] 
%     %...............................................................................
%     d1uxdx = uxj .* (  -(xj - sstj)/gammay^2);
%     d2uxdx = uxj .* ( (-(xj - sstj)/gammay^2).^2 - (1/gammay^2) );
%     % Third and fourth derivatives added, yj transformed to xj (Le Gland, 25/04/2019)
%     d3uxdx = uxj .* ( (-(xj - sstj)/gammay^2).^3 + (xj - sstj)/gammay^4 + 2*(xj - sstj)/gammay^4 );
%     d4uxdx = uxj .* ( (-(xj - sstj)/gammay^2).^4 - 3*((xj - sstj).^2)/gammay^6 - 3*((xj - sstj).^2)/gammay^6 + (3/gammay^4) );
%     %...............................................................................
%     %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% % $$$     %...............................................................................
% % $$$     qyj = q10; %Control case without gaussian. 
% % $$$     uxj = mup .* qyj .* (din ./ (din + knp0)); %[d-1] 
% % $$$     %...............................................................................
% % $$$     d1uxdx  = 0d0;
% % $$$     d2uxdx = 0d0;
% % $$$     sigmax = 0d0;
% % $$$     %...............................................................................
%     %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%     %===============================================================================
% end

% Scheme to "artificially" increase competition (Le Gland, 14/05/2019)
% d1uxdx =  5*d1uxdx;
% d1gxdx =  5*d1gxdx;
% d2uxdx = 25*d2uxdx;
% d2gxdx = 25*d2gxdx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ORDINARY DIFFERENTIAL EQUATIONS (ODEs):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
% $$$ ux = uxj + (1/2)*(sigmaxj.^2).*(d2uxdx); %Community upake rate gain.
%...................................................................................
% ux = uxj + (1/2)*(sigmaxj.^2).*(d2uxdx + numut.*d4uxdx) - (1/2)*3*numut.*d2uxdx; %Community upake rate gain [d-1].
%SIMPLIFICATION (Le Gland, 08/01/2019)
% ux = uxj + (1/2)*(sigmaxj.^2).*d2uxdx;
% ux = uxj + (1/2)*(sigmaxj.^2).*(d2uxdx + numut.*d4uxdx) - (1/2)*3*numut.*d2uxdx + (3/24)*(sigmaxj.^4).*d4uxdx;
% Add kurtosis term (Le Gland, 25/04/2019)
% ux = ux + (3/24)*(sigmaxj.^4).*(d4uxdx - d4gxdx);
% Case with 2 traits (Le Gland, 03/07/2019)
%uxy = uxyj + (1/2)*sigma(i)*d2uxy ; % Need matricial writing
uxy = uxyj + (1/2)*(xxvar.*d2uxydxdx + yyvar.*d2uxydydy) + xycov.*d2uxydxdy;

%Take exudation into account (Le Gland, 03/10/2019)

if strcmp(keyAssimConstant,'not')

    betap = betap_max .* exp( -1/2 .* (xmj - betap_xmax).^2 ./ (betap_xrange.^2) );
    d1betapdx   = betap .* ( - (xmj - betap_xmax) ) ./ (betap_xrange.^2);
    d2betapdxdx = betap .* ( (xmj - betap_xmax).^2 ./ (betap_xrange.^4) - 1 ./ (betap_xrange.^2) );

    exyj      = (1-betap).*uxyj;
    d1exydx   = (1-betap).*d1uxydx - d1betapdx.*uxyj;
    % d1exydy   = 0;
    d2exydxdx = (1-betap).*d2uxydxdx - 2*d1betapdx.*d1uxydx - d2betapdxdx.*uxyj;
    % d2exydydy = 0;
    % d2exydxdy = 0;
    % exy = exyj + (1/2)*(xxvar.*d2exydxdx + yyvar.*d2exydydy) + xycov.*d2exydxdy;
    exy = exyj + (1/2)*(xxvar.*d2exydxdx); % Is there a first derivative term ? Probably not

else
 
    betap = betap_max;
    %d1exydx = 0;
    %d2exydxdx = 0;
    % ERROR: d1exydx and d2exydxdx have to be proportional to the
    % derivatives of uxyj (Le Gland, 04/11/2019)
    d1exydx   = (1-betap).*d1uxydx;
    d2exydxdx = (1-betap).*d2uxydxdx;
    exy = (1-betap)*uxy;
    
end

% Exudation derivatives for the remaining terms (Le Gland, 04/11/2019)
d1exydy   = (1-betap).*d1uxydy;
d2exydydy = (1-betap).*d2uxydydy;
d2exydxdy = (1-betap).*d2uxydxdy;

%...................................................................................
%===================================================================================
%...................................................................................
% $$$ gx = gxj + (1/2)*(sigmaxj.^2).*(d2gxdx); %Community grazing rate loss [d-1]. %WRONG!!!!
%...................................................................................
%gx = (1./phy) .* Qfeeding .* Vmax; %Community grazing rate loss [d-1] %OKAY (NOTE: Ptot = phy)
%gxy = (1./phy) .* Qfeeding .* Vmax; % 2 traits (Le Gland, 16/07/2019)
%Correction term for high grazing rate on small phytoplankton (Le Gland, 26/08/2020)
gxy = (1./phy) .* Qfeeding .* Vmax .* exp(amuz^2 .* invsigma(:,2,2) ./ (2*galfa));

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
%Fphy = ux .* phy; %Phy primary production [mmolN * m-3 * d-1]
Fphy = uxy .* phy; % 2 traits (Le Gland, 16/07/2019)
%Gphy = gx .* phy; %Phy grazing mortality [mmolN * m-3 * d-1]
Gphy = gxy .* phy; % 2 traits (Le Gland, 16/07/2019) 
% Mphy = mp  * phy; %Phy natural mortality [mmolN * m-3 * d-1] %Linear.
Mphy = mp * phy .* (Q10h.^((tempj - temp0)/10)); % Temperature dependence on mortality (Le Gland, 11/11/2019)
% Mphy = mp  * (phy.*phy); %Phy natural mortality [mmolN * m-3 * d-1] %Quadratic.
% Mphy = mp  * (phy.*phy) .* (Q10h.^((tempj - temp0)/10));
% Test with grazing equal to growth (Le Gland, 09/05/2019)
% Mphy = 0;
%...................................................................................
Fzoo = Gphy; %Zoo second production [mmolN * m-3 * d-1]
Ezoo = (1-betaz) * Fzoo; %Zoo exudation [mmolN * m-3 * d-1] 
%Ephy = (1-betap) * Fphy; %Phy exudation [mmolN * m-3 * d-1]
Ephy = exy .* phy; % Variable exudation rate (Le Gland, 03/10/2019)
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
%.......................................  ............................................
%ORIGINAL:
% dXAVEdt = (sigmaxj.^2) .* (d1uxdx - d1gxdx + numut.*d3uxdx) - 3*numut.*d1uxdx; 
% dXVARdt = (sigmaxj.^4) .* (d2uxdx - d2gxdx + numut.*d4uxdx) - (sigmaxj.^2) .* (5*numut.*d2uxdx) + (2*numut.*uxj); 
%...................................................................................
%SIMPLIFIED: DOES NOT WORK!!! BLOWS UP AT TIME 700 days.
% dXAVEdt = (sigmaxj.^2) .* (d1uxdx - d1gxdx);
% dXVARdt = (sigmaxj.^4) .* (d2uxdx - d2gxdx) - (sigmaxj.^2) + (2*numut.*uxj); 
%...................................................................................
%SECOND SIMPLIFICATION (Le Gland, 08/01/2019)
% dXAVEdt = 0; % No variance case
% dXAVEdt = (sigmaxj.^2) .* (d1uxdx - d1gxdx);
% dXAVEdt = (sigmaxj.^2) .* (d1uxdx - d1gxdx + numut.*d3uxdx) - 3*numut.*d1uxdx + (1/2)*(sigmaxj.^4).*d3uxdx;
% Add kurtosis term (Le Gland, 25/04/2019)
% dXAVEdt = dXAVEdt + (3/6)*(sigmaxj.^4) .* (d3uxdx - d3gxdx);
% dXVARdt = 0; % Constant variance case
% dXVARdt = (sigmaxj.^4) .* (d2uxdx - d2gxdx) + (sigmaxj.^2) .* (numut.*d2uxdx) + (2*numut.*uxj);
% dXVARdt = (sigmaxj.^4) .* (d2uxdx - d2gxdx + numut.*d4uxdx) - (sigmaxj.^2) .* (5*numut.*d2uxdx) + (2*numut.*uxj);
%...................................................................................
% Case with 2 traits (Le Gland, 03/07/2019)
% dXAVEdt = V * (d1uxy - d1gxy);
% dXVARdt = V * (d2uxy - d2gxy) * V + 2 * Mut * (uxj + 1/2*(sum(sum(V.*d2uxy))));

% Intent to take betap into account (Le Gland, 20/09/2019)
% uxyj    = betap * uxyj;
% d1uxydx = betap * d1uxydx;
% d1uxydy = betap * d1uxydy;
% d2uxydxdx = betap * d2uxydxdx;
% d2uxydydy = betap * d2uxydydy;
% d2uxydxdy = betap * d2uxydxdy;
% d1gxydx = 0.4*d1gxydx;
% d1gxydy = 0.4*d1gxydy;
% d2gxydxdx = 0.4*d2gxydxdx;
% d2gxydydy = 0.4*d2gxydydy;
% d2gxydxdy = 0.4*d2gxydxdy;
%xycov = 0;

% Non matricial writing (Le Gland, 12/07/2019)
% Exudation added on 03/10/2019
dXAVEdt = xxvar .* (d1uxydx - d1gxydx - d1exydx) + xycov .* (d1uxydy - d1gxydy - d1exydy);
dYAVEdt = xycov .* (d1uxydx - d1gxydx - d1exydx) + yyvar .* (d1uxydy - d1gxydy - d1exydy);

% Replaced numut by numutx and numuty (Le Gland, 10/09/2019)
% Exudation added on 03/10/2019
dXXVARdt = (xxvar.^2) .* (d2uxydxdx - d2gxydxdx - d2exydxdx) + (xxvar.*xycov) .* (d2uxydxdy - d2gxydxdy - d2exydxdy) + ...
           (xycov.^2) .* (d2uxydydy - d2gxydydy - d2exydydy) + numutx .* (2.*(uxy-exy)); % (2.*uxyj + xxvar.*d2uxydxdx + ...
           %yyvar.*d2uxydydy + 2.*xycov.*d2uxydxdy);
dYYVARdt = (xycov.^2) .* (d2uxydxdx - d2gxydxdx - d2exydxdx) + (yyvar.*xycov) .* (d2uxydxdy - d2gxydxdy - d2exydxdy) + ...
           (yyvar.^2) .* (d2uxydydy - d2gxydydy - d2exydydy) + numuty .* (2.*(uxy-exy)); % (2.*uxyj + xxvar.*d2uxydxdx + ...
           %yyvar.*d2uxydydy + 2.*xycov.*d2uxydxdy);
dXYCOVdt = (xxvar.*xycov) .* (d2uxydxdx - d2gxydxdx - d2exydxdx) + (yyvar.*xycov) .* (d2uxydydy - d2gxydydy - d2exydydy) + ...
           (xxvar.*yyvar + xycov.^2) .* (d2uxydxdy - d2gxydxdy - d2exydxdy) + dXYCOVdt_control;
% dXYCOVdt = 0;

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
dXAVE_STARdt = (dPHYdt .* xave) + (dXAVEdt .* phy);
% dXVAR_STARdt = (dPHYdt .* xvar) + (dXVARdt .* phy);
% Change by Le Gland (23/04/2019)
dXXVAR_STARdt = (dPHYdt .* xxvar) + (dXXVARdt .* phy) + (dPHYdt .* xave .* xave) + 2*(dXAVEdt .* xave .* phy);
% dXSTD_STARdt = (dPHYdt .* xstd) + (dXSTDdt .* phy);
% Case with 2 traits (Le Gland, 03/07/2019)
dYAVE_STARdt  =  (dPHYdt .* yave)  + (dYAVEdt .* phy);
dYYVAR_STARdt =  (dPHYdt .* yyvar) + (dYYVARdt .* phy) + (dPHYdt .* yave .* yave) + 2*(dYAVEdt .* yave .* phy);
dXYCOV_STARdt =  (dPHYdt .* xycov) + (dXYCOVdt .* phy) + (dPHYdt .* xave .* yave) + (dXAVEdt .* yave .* phy) + (dYAVEdt .* xave .* phy);
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
%dPONdt(1)
%DIFFpon(1)
%ADVpon(1)
%...................................................................................
boxdot = dBOXdt; 
%...................................................................................
%===================================================================================
if strcmp(keyPhysics,'not')
    %...............................................................................
    xavedot  = dXAVEdt;
    xxvardot = dXXVARdt;
    % xstddot = dXSTDdt;
    % Case with 2 traits (Le Gland, 03/07/2019)
    yavedot  = dYAVEdt;
    yyvardot = dYYVARdt;
    xycovdot = dXYCOVdt;
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
    xave_stardot  = dXAVE_STARdt  + DIFFxave_star;
    xxvar_stardot = dXXVAR_STARdt + DIFFxxvar_star;
    % xstd_stardot = dXSTD_STARdt + DIFFxstd_star;
    % Case with 2 traits (Le Gland, 03/07/2019)
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
    %Vdot = [xavedot;xvardot;xstddot;phydot;zoodot;dindot;pondot;boxdot];
    % 2 traits (Le Gland, 12/07/2019)
    % Vdot = [xavedot;yavedot;xxvardot;yyvardot;xycovdot;phydot;zoodot;dindot;pondot;boxdot];
    Vdot = [phydot;zoodot;dindot;pondot;boxdot;xavedot;yavedot;xxvardot;yyvardot;xycovdot];
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    %Vdot = [xave_stardot;xvar_stardot;xstd_stardot;phydot;zoodot;dindot;pondot;boxdot];
    % 2 traits (Le Gland, 12/07/2019)
    % Vdot = [xave_stardot;yave_stardot;xxvar_stardot;yyvar_stardot;xycov_stardot;phydot;zoodot;dindot;pondot;boxdot];
    Vdot = [phydot;zoodot;dindot;pondot;boxdot;xave_stardot;yave_stardot;xxvar_stardot;yyvar_stardot;xycov_stardot];
    %...............................................................................
end

return
%...................................................................................
%===================================================================================
%***********************************************************************************
