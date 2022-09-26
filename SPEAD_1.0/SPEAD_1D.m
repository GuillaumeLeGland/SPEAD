%********************************************************************
% SPEAD (Simulating Plankton Evolution with Adaptive Dynamics) version 1.0
% Code written by Guillaume Le Gland and Sergio M. Vallina
% Published on the 20th of October 2020
% Presented in Le Gland et al., 2020, under review in GMD.
%********************************************************************

%********************************************************************
% PROGRAM: SPEAD_1D.M
%********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FILE HEADER -- LOAD PACKAGES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************************************************************
% more off
% close all
% clear all
% format short g
%--------------------------------------------------------------------
% addpath(genpath('TOOLBOX/'));
%--------------------------------------------------------------------
%====================================================================
%--------------------------------------------------------------------

function [dat14,dat15,dat2020a,dat1010,dat1020,dat1022,dat1030,...
    dat2010,dat2020,dat2030,dat70,dat80,dat24]=SPEAD_1D()

[mypackages] = myheadloadpackages(false); %Structure array to pass on my head pkg as input argument to functions.
%--------------------------------------------------------------------
%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal;
%--------------------------------------------------------------------
%====================================================================
%********************************************************************
global galfa gbeta
global gzmax kgz mz betaz mpower 
global mp Isat InhFac numutx numuty
global alp0 mup0
global aalp amup aknp
global Q10a Q10h % Distinct partition coefficients for auto and heterotrophic processes
global temp0 temp
global tcounter jcounter
global keyPhysics keySinking
global omePhy epsZoo omeZoo md % epsPhy (not required because we neglect phytoplankton exudation) 
global t0 deltat ndays nyear tmax tspan 
global zdepths ndepths deltaz
global Ixave Iyave Ixxvar Iyyvar Ixycov
global Ixave_K Ixxvar_K Iyave_T Iyyvar_T
global Iphy Izoo Idin Ipon Ibox %continuous model
global Jphy Jzoo Jdin Jpon Jbox %discrete model  
global nxphy nyphy nzoo ndin npon nbox % Number of x-trait and y-trait values in the discrete model + number of zoo, din and pon
global mtot0 %discrete model
global KZ
global parz0
global kw wsink
global jday
global keyTraitDiffusion
global xaxis yaxis xdel ydel 
%...................................................................................
global UXYout GXYout % 2-trait model outputs
global UXout GXout UYout GYout % 1-trait model outputs
%...................................................................................
global uphydaydisc gphydaydisc %OUTPUTS FROM ODE45 OF DISCRETE MODEL.  
%...................................................................................
global FPHYdaydisc GPHYdaydisc %OUTPUTS FROM ODE45 OF DISCRETE MODEL. 
%...................................................................................

tic 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL KEYS AND PARAMETERS PACKAGES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%KEYS
%........................................................................
[keyKN,keyTOPT,key2T,keyDisc,keyModelResol,keyPhysics,keySinking,keyPARseasonality,...
keyPARextinction,keyKTW,keyTraitDiffusion] = ...
SPEAD_1D_keys;
%........................................................................
%========================================================================
%PARAMETERS
%........................................................................
[deltat,nyear,ndays,deltaz,ndepths,betaz,gzmax,kgz,mz,mpower,Q10a,Q10h,temp0,Isat,InhFac,mup0,alp0,amup,...
aalp,aknp,mp,epsZoo,omePhy,omeZoo,md,galfa0,gbeta0,kw,wsink,xmin,xmax,xsigmaPcnt,ymin,ymax,ysigmaPcnt,...
nxphy,nyphy,numutX0,numutY0,phy0,zoo0,din0,pon0,box0] = ...
SPEAD_1D_parameters(keyKTW,keyTraitDiffusion,keyModelResol);
%........................................................................

%%%%%%%%%%%%%%%%%%%%%
%TEMPORAL RESOLUTION:
%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
t0 = deltat;
tmax = ndays*nyear; 
tspan = t0:deltat:tmax; 
%........................................................................
nsteps = length(tspan); % Should be equal to (ndays/deltat)*nyear
ttmax = ndays; %[days]
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%
%VERTICAL RESOLUTION:
%%%%%%%%%%%%%%%%%%%%%
%========================================================================

% Change in the depth levels
% First level begins at deltaz/2 (box from0 to deltaz)
zmin = 0;
zmax = ndepths*deltaz;
zdepths = zmin+deltaz/2:deltaz:zmax-deltaz/2;

%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%
%EXTERNAL FORCINGS:
%%%%%%%%%%%%%%%%%%%

[itemp,iparz0,PAR2D,imld,iKZ] = SPEAD_1D_externalForcings(ndepths,zdepths,ndays,kw);

[CHL_obs,PP_obs,NO3_obs,PON_obs] = SPEAD_1D_obsLoad();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAKE LONGER ARRAYS (mld, par, KZ, KZI) ACCORDING TO NUMBER OF YEARS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%........................................................................
if ndepths > 1
    KZ  = repmat(iKZ,[1,nyear]);
end
mld = repmat(imld,[1,nyear]);
temp = repmat(itemp,[1,nyear]);
parz0 = repmat(iparz0,[1,nyear]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STABILITY CONDITION FOR TURBULENCE SCHEME: 
%(not really relevant when using "ode45" solver 
% because it uses a variable adaptive dt -- 
% use then "ode4" which has aconstant dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------
%ONLY NEEDED WHEN USING EXPLICIT SCHEMES -- 
% NOT NEEDED WHEN USING IMPLICIT SCHEMES
%======================================
%--------------------------------------
% dtdiff < min(dz^2/(2*kzI))
%--------------------------------------
%........................................................................
% Stability check is not necessary if the vertical diffusion scheme is
% implicit or has a split time step
% if ndepths > 1
%     maxKZ = max(KZ);
%     stability = min(deltaz^2./(2*maxKZ)) % Maximum dt allowed.
%     %........................................................................
%     if (deltat*2) >= stability 
%         disp('Be careful: Condition of stability violated!')
%         disp('reduce dt or increase dz')
%         pause
%     end
%     %........................................................................
%     dtadv = min((deltaz^2)./((wsink*deltaz) + (2*maxKZ))); %max time step from von Neumann stability analysis
%     %........................................................................
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHY HALF-SAT RESOLUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
% Define Knmin and Knmax based on the x trait
Knmin = exp(xmin);
Knmax = exp(xmax);
xdel = (xmax-xmin)/(nxphy-1); % Probably the most consistent expression, since no point has to be removed (Le Gland, 04/06/2019)   
%...................................................................................
xaxis  = xmin:xdel:xmax; %[log(um)] 
%...................................................................................
Kn = exp(xaxis); %Equivalent Spherical Diameter [um]
%...................................................................................
show_Kn = [[1:nxphy]',xaxis(:),Kn(:)]
%...................................................................................
ydel = (ymax-ymin)/(nyphy-1);
yaxis = ymin:ydel:ymax;
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IF YOU WANT TO REMOVE ALL PHYSICAL PROCESSES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
if strcmp(keyPhysics,'not') && ndepths > 1
    KZ (:,:) = 0; %no turbulence
end
%........................................................................
if strcmp(keySinking,'not')
    wsink = 0; %no sinking
end
%........................................................................
if strcmp(keyPARextinction,'not')
    kw = 0; %no light extinction in water column.
end
%........................................................................
if strcmp(keyPARseasonality,'not')
   iparz0 = Isat*ones(1,ndays); %Constant PAR at the surface.
    parz0 = Isat*ones(1,ndays*nyear); %Constant PAR at the surface.
    PAR2D = exp(-kw*zdepths(:))*iparz0;
end

dat14={itemp,iparz0,PAR2D,imld,iKZ,14,mypackages};
dat15={mup0,amup,[0.1,0.5,2.0],temp0,Q10a,18:4:30,15};

%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIAL CONDITIONS FOR THE CONTINOUS TRAIT MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...............................................................................
xave0  = (xmax + xmin) / 2; %Phy diameter average [log(um)]
xstd0  = xsigmaPcnt * (xmax - xmin);
xxvar0 = xstd0^2; % Le Gland, 04/07/2019
%...............................................................................
yave0 = (ymax + ymin) / 2;
ystd0 = ysigmaPcnt * (ymax - ymin);
yyvar0 = ystd0^2;
% Correlation between traits is initialized at zero 
xycov0 = 0;
%...............................................................................
%===================================================================================
% The raw moments multiplied by phy0 ("star") are the conservative
% quantities to be used in the vertical mixing function
if strcmp(keyPhysics,'yes')
    %...............................................................................
    xave_star0  = xave0*phy0;
    xxvar_star0 = (xxvar0+xave0.^2)*phy0;
    %...............................................................................
    yave_star0  = yave0*phy0;
    yyvar_star0 = (yyvar0+yave0.^2)*phy0; 
    %...............................................................................
    xycov_star0 = (xave0.*yave0)*phy0;
    %...............................................................................
end

%===================================================================================
%...................................................................................
%CONSTANT VERTICAL PROFILES OF THE INITIAL CONDITIONS:
%........................................................................
vdepths=ones(ndepths,1); %column vector of ones with ndepths vertical nodes.
%........................................................................
% Add a second trait (Le Gland, 01/07/2019)
if strcmp(keyPhysics,'not') || ndepths == 1
    %....................................................................
    PDFxave0  = xave0 *vdepths;
    PDFxxvar0 = xxvar0*vdepths; % Le Gland, 04/07/2019
    PDFyave0  = yave0 *vdepths;
    PDFyyvar0 = yyvar0*vdepths;
    PDFxycov0 = xycov0*vdepths;
    %....................................................................
elseif strcmp(keyPhysics,'yes')
    %....................................................................
    PDFxave_star0  = xave_star0 *vdepths;
    PDFxxvar_star0 = xxvar_star0*vdepths; % Le Gland, 04/07/2019
    PDFyave_star0  = yave_star0 *vdepths;
    PDFyyvar_star0 = yyvar_star0*vdepths;
    PDFxycov_star0 = xycov_star0*vdepths;
    %....................................................................
end
%........................................................................
PHY0 = phy0*vdepths; %column vector with I.C. for Phy (ndepths vertical nodes).
ZOO0 = zoo0*vdepths; %column vector with I.C. for Zoo (ndepths vertical nodes).
DIN0 = din0*vdepths; %column vector with I.C. for Nut (ndepths vertical nodes).
PON0 = pon0*vdepths; %column vector with I.C. for Det (ndepths vertical nodes).
BOX0 = box0*vdepths; %column vector to check mass balance
%........................................................................
gbeta = gbeta0;
galfa = galfa0*vdepths; %Constant grazing-alfa KTW parameter.
numutx = numutX0*vdepths; %column vector for x mutation rate
numuty = numutY0*vdepths; %column vector for y mutation rate
%........................................................................
%========================================================================
%........................................................................
%COLUMN VECTOR OF INITIAL CONDITIONS:
%........................................................................
if strcmp(keyPhysics,'not')
    %....................................................................
    PDF = [PDFxave0,PDFyave0,PDFxxvar0,PDFyyvar0,PDFxycov0]; % 2-trait aggregate model
    PDF_K = [PDFxave0,PDFxxvar0]; % 1-trait aggregate model with Kn as a trait
    PDF_T = [PDFyave0,PDFyyvar0]; % 1-trait agggregate model with Topt as a trait
    %....................................................................
elseif strcmp(keyPhysics,'yes')
    %....................................................................
    PDF = [PDFxave_star0,PDFyave_star0,PDFxxvar_star0,PDFyyvar_star0,PDFxycov_star0]; % 2-trait aggregate model
    PDF_K = [PDFxave_star0,PDFxxvar_star0]; % 1-trait aggregate model with Kn as a trait
    PDF_T = [PDFyave_star0,PDFyyvar_star0]; % 1-trait agggregate model with Topt as a trait
    %....................................................................
end
%........................................................................
BIO = [PHY0,ZOO0,DIN0,PON0,BOX0]; %Biomasses of plankton and nutrients.
%........................................................................
PDF_V0 = PDF(:); % Transform a matrix (ndepths, 5) into a vector (ndepths * 5)
PDF_V0_K = PDF_K(:);
PDF_V0_T = PDF_T(:);
%........................................................................
BIO_V0 = BIO(:);
%........................................................................
V0 = [BIO_V0;PDF_V0];
V0_K = [BIO_V0;PDF_V0_K];
V0_T = [BIO_V0;PDF_V0_T];
%........................................................................
ntot0 = sum(BIO_V0(:)); %initial total mass over the column water [mmolN*m-3]
%........................................................................
ntot0pernode = ntot0/ndepths; %initial total mass at each grid cell [mmolN*m-3]
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIAL CONDITIONS FOR THE DISCRETE TRAIT MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------------
%........................................................................
nzoo =  1; %Number of zoo species (must be always only one!)
ndin =  1; %DIN (must be always only one!)
npon =  1; %PON (must be always only one!)
nbox =  1; %Virtual box (must be always only one!)
%........................................................................
%========================================================================
%DISCRETE TRAITS:
%........................................................................
xsigma = xsigmaPcnt*(xmax-xmin);
xmean = (xmin + xmax) / 2;
xplus = xmean + xsigma;
xless = xmean - xsigma;
Knplus = exp(xplus);
Knless = exp(xless);
%........................................................................
ymean = (ymin + ymax) / 2;
ysigma = ysigmaPcnt*(ymax-ymin); 
yplus = ymean + ysigma;
yless = ymean - ysigma;
%........................................................................
xdelta = xplus - xless; 
ydelta = yplus - yless;
%........................................................................
%========================================================================
%........................................................................
fxtraitphy = (1.0 / (xsigma * sqrt(2*pi))) * exp( -(xaxis - xmean).^2 / (2*xsigma^2) );
fytraitphy = (1.0 / (ysigma * sqrt(2*pi))) * exp( -(yaxis - ymean).^2 / (2*ysigma^2) );
% fxytraitphy is a multivariate normal distribution without correlation (Le Gland, 04/07/2019)    
fxytraitphy = fxtraitphy'*fytraitphy;
%........................................................................
%========================================================================

dat2020a={xaxis,fxtraitphy,xless,xplus,xmin,xmax,Knless,Knplus,Knmin,Knmax};

%========================================================================
%FOR GAUSSIAN DISTRIBUTION OF INITIAL CONDITIONS OF PHYTOPLANKTON: 
%........................................................................
phydisc0gaussian = max(eps, phy0 * (fxytraitphy .* (xdel*ydel))); % "eps" is to avoid strictly zero values
PHYdisc0gaussian1D = zeros(ndepths,nxphy,nyphy);
for i = 1:ndepths 
    PHYdisc0gaussian1D(i,:,:) = phydisc0gaussian;
end
PHYdisc0 = PHYdisc0gaussian1D(:); % Turns the trait distribution into a vector
%........................................................................
ZOOdisc0 = zoo0 * ones(ndepths*nzoo,1);
DINdisc0 = din0 * ones(ndepths*ndin,1);
PONdisc0 = pon0 * ones(ndepths*npon,1);
BOXdisc0 = box0 * ones(ndepths*nbox,1);
%........................................................................
%========================================================================
%COLUMN VECTOR OF INITIAL CONDITIONS:
%........................................................................
Vdisc0 = [PHYdisc0;ZOOdisc0;DINdisc0;PONdisc0;BOXdisc0];
%........................................................................
mtot0 = sum(Vdisc0);
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE THE ROWS CORRESPONDING TO EACH STATE-VARIABLE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%FOR CONTINOUS TRAIT MODEL:
%------------------------------------------------------------------------
%[NOTE: Iphy, Izoo, Idin, Ipom]
%------------------------------------------------------------------------
%........................................................................
mlength = length(BIO_V0);
%........................................................................
% These indices are common to the 1-trait and 2-trait models
Iphy = 1:ndepths;
Izoo = (ndepths*1)+1:(ndepths*2);
Idin = (ndepths*2)+1:(ndepths*3);
Ipon = (ndepths*3)+1:(ndepths*4);
Ibox = (ndepths*4)+1:(ndepths*5);
%........................................................................
% Indices for the 2-trait model
Ixave  = mlength + (1:ndepths);
Iyave  = mlength + ((ndepths*1)+1:(ndepths*2));
Ixxvar = mlength + ((ndepths*2)+1:(ndepths*3));
Iyyvar = mlength + ((ndepths*3)+1:(ndepths*4));
Ixycov = mlength + ((ndepths*4)+1:(ndepths*5));
% Indices in 1-trait models (Le Gland, 25/11/2019)
Ixave_K  = mlength + (1:ndepths);
Iyave_T  = mlength + (1:ndepths);
Ixxvar_K = mlength + ((ndepths*1)+1:(ndepths*2));
Iyyvar_T = mlength + ((ndepths*1)+1:(ndepths*2));
%........................................................................
wIndex = [Iphy;Izoo;Idin;Ipon;Ibox;Ixave;Iyave;Ixxvar;Iyyvar;Ixycov]
%........................................................................
%========================================================================
%FOR DISCRETE TRAIT MODEL:
%------------------------------------------------------------------------
%[NOTE: Jphy, Jzoo, Jdin, Jpom]
%------------------------------------------------------------------------
%........................................................................
%nvar = [0,nphy,nzoo,ndin,npon,nbox];
nvar = [0,nxphy*nyphy,nzoo,ndin,npon,nbox]; % Case with 2 traits (Le Gland, 05/07/2019)
%........................................................................
varname{01} = 'phy';
varname{02} = 'zoo';
varname{03} = 'din';
varname{04} = 'pon';
varname{05} = 'box';
%........................................................................
nsum1 = 0; 
nsum2 = 0; 
%........................................................................
for j = 1:length(nvar)-1
    %....................................................................
    nsum1 = nsum1 + nvar(j);
    nsum2 = nsum2 + nvar(j+1);
    %....................................................................
    Jvarj = (ndepths*(nsum1))+1:(ndepths*(nsum2));
    %....................................................................
    Jstr = ['J',varname{j}];
    %....................................................................
    %% assignin('base',Jstr,Jvarj); %DO NOT USE THIS OLD APPROACH
    %....................................................................
    myassign = [Jstr,' = Jvarj;']; %BETTER TO USE THIS NEW APPROACH FOR GLOBAL VARIABLES
    %....................................................................
    eval(myassign)
    %....................................................................
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRAIT DIFFUSION FOR DISCRETE MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
if strcmp(keyTraitDiffusion,'yes')
    mup_mean = mup0 .* exp(amup*mean(xaxis)) .* Q10a^(mean(yaxis)/10); % [d-1]
    nux = numutX0;
    nuy = numutY0;
    kppx = nux.*mup_mean;
    kppy = nuy.*mup_mean;
    maxKPPx = max(kppx(:));
    maxKPPy = max(kppy(:));
    %...............................................................................
    %===============================================================================
    %...............................................................................
    StabilityKPP  = min(xdel^2./(2*maxKPPx),ydel^2./(2*maxKPPy))
    %...............................................................................
    if deltat < StabilityKPP

        checkStabilityDIFF = 'StableNumerically ==> OK';
        checkStabilityDIFF

    elseif deltat >= StabilityKPP

        checkStabilityDIFF = 'UnstableNumerically ==> ERROR!!!';
        checkStabilityDIFF

        wStability = [deltat,StabilityKPP] 

        disp('Error: Violation of numerical stability for trait diffusion!!')
        disp('To avoid this problem, either decrease "dt" or increase "dx"!!')
        pause
    end
end
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RUNGE-KUTTA ODE45 SOLVER FOR CONTINOUS TRAIT MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
display('*** continuous model ode45 solving -- START ***')
%...................................................................................
jday  = 1;
%...................................................................................
tcounter = 0;
jcounter = 0;
%...................................................................................
Voutcont = zeros(length(tspan),length(V0));
%...................................................................................
ode45options = odeset('AbsTol',1e-12,'RelTol',1e-6);
%...................................................................................
toc

% Aggregate model 2 traits
if strcmp(key2T,'yes')
    tic
    [Voutcont] = ode4(@SPEAD_gaussecomodel1D_ode45eqs,tspan,V0);
     Voutcont = Voutcont';
    toc
end
%...................................................................................
% Outputs from 1-trait models
if strcmp(keyKN,'yes')
    tic
    [Voutcont_K] = ode4(@SPEAD_gaussecomodel1D_ode45eqs_Ksat,tspan,V0_K);
     Voutcont_K = Voutcont_K';
    toc
end
if strcmp(keyTOPT,'yes')
    tic
    [Voutcont_T] = ode4(@SPEAD_gaussecomodel1D_ode45eqs_Topt,tspan,V0_T);
     Voutcont_T = Voutcont_T';
    toc
end
%...................................................................................
% Option to save all variables
%% save('SPEAD_1D_Voutcont.mat','Voutcont')
%...................................................................................
display('*** continuous model ode45 solving -- END ***')
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RUNGE-KUTTA ODE45 SOLVER FOR DISCRETE TRAIT MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
display('*** discrete model ode45 solving -- START ***')
%...................................................................................
jday  = 1;
%...................................................................................
icounter = 1;
%...................................................................................
Voutdisc = zeros(length(tspan),length(Vdisc0));
%...................................................................................
ode45options = odeset('AbsTol',1e-12,'RelTol',1e-6);
%...................................................................................
if strcmp(keyDisc,'yes')
    tic
    [Voutdisc] = ode4(@SPEAD_discretemodel1D_ode45eqs,tspan,Vdisc0);
    toc
    Voutdisc = Voutdisc';
end
%...................................................................................
%% save('SPEAD_1D_Voutdisc.mat','Voutdisc')
%...................................................................................
display('*** discrete model ode45 solving -- END ***')
%...................................................................................
toc
%===================================================================================
%***********************************************************************************
pause(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STATE VARIABLES CONTINOUS TRAIT MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%GET DAILY AVERAGES:
%...................................................................................
if strcmp(key2T,'yes')
    [Vodecont,todecont] = SPEAD_1D_dailyAve(Voutcont,deltat);
end
if strcmp(keyKN,'yes')
    [Vodecont_K,todecont] = SPEAD_1D_dailyAve(Voutcont_K,deltat);
end
if strcmp(keyTOPT,'yes')
    [Vodecont_T,todecont] = SPEAD_1D_dailyAve(Voutcont_T,deltat);
end
%...................................................................................
SPEAD_1D_analysis % Script to analyze the model outputs
%...................................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT STATISTICS OF BOTH CONTINUOUS AND DISCRETE MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%CONTINOUS:
if strcmp(key2T,'yes')
    %...................................................................................
    fignum = 1010;
    dat1010 = {MUPsspcont,MUZsspcont,NTOTsspcont,ndepths,ndays,MUPmin,MUPmax,...
        MUZmin,MUZmax,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages};
    %...................................................................................
    fignum = 1020;
    dat1020 = {NTOTsspcont,CHLsspcont,PPsspcont,ZOOsspcont,DINsspcont,PONsspcont,ndepths,ndays,PHYmin,PHYmax,...
        ZOOmin,ZOOmax,DINmin,DINmax,PONmin,PONmax,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages};
    %...................................................................................
    % Compare model and observations
    fignum = 1022;
    dat1022 = {12*(106/16)*PP_obs,CHL_obs,NO3_obs,PON_obs,12*(106/16)*PPsspcont,...
        CHLsspcont,DINsspcont,PONsspcont,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages};
    %...................................................................................
    fignum = 1030;
    dat1030 = {logESDphysspAveCont,logESDphysspStdCont,TOPTphysspAveCont,TOPTphysspStdCont,physspCorCont,PHYTsspcont,...
        ndepths,ndays,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,logESDaveMax,logESDaveMin,logESDstdMax,logESDstdMin,...
        TOPTaveMax,TOPTaveMin,TOPTstdMax,TOPTstdMin,CorrelationAbsMax,fignum,mypackages};
    %...................................................................................
end
%===================================================================================
%DISCRETE:
if strcmp(keyDisc,'yes')
    %...................................................................................
    fignum = 2010;
    dat2010 = {MUPsspdisc,MUZsspdisc,NTOTsspdisc,ndepths,ndays,MUPmin,MUPmax,...
        MUZmin,MUZmax,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages};
    %...................................................................................
    fignum = 2020;
    dat2020 = {NTOTsspdisc,CHLsspdisc,PPsspdisc,ZOOsspdisc,DINsspdisc,PONsspdisc,ndepths,ndays,PHYmin,PHYmax,...
        ZOOmin,ZOOmax,DINmin,DINmax,PONmin,PONmax,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages};
    %...................................................................................
    fignum = 2030;
    dat2030 = {logESDphysspAveDisc,logESDphysspStdDisc,TOPTphysspAveDisc,TOPTphysspStdDisc,physspCorDisc,PHYTsspdisc,...
        ndepths,ndays,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,logESDaveMax,logESDaveMin,logESDstdMax,logESDstdMin,...
        TOPTaveMax,TOPTaveMin,TOPTstdMax,TOPTstdMin,CorrelationAbsMax,fignum,mypackages};
    %...................................................................................
end
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTS OF TRAIT DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
if strcmp(key2T,'yes') && strcmp(keyDisc,'yes') && strcmp(keyModelResol,'1D')
    fignum = [70];
    dat70 = {PHYTsspcont,logESDphysspAveCont,logESDphysspStdCont,TOPTphysspAveCont,TOPTphysspStdCont,physspCorCont,PHYsspdisc3D,xaxis,yaxis,itemp,DINsspdisc,fignum};
end
%...................................................................................
if strcmp(key2T,'yes') && strcmp(keyKN,'yes') && strcmp(keyTOPT,'yes')
    fignum = [80];
    dat80 = {DINsspcont,DINsspcont_K,temp(:,1:360),...
    PHYTsspcont,PHYTsspcont_K,PHYTsspcont_T,XAVE_sspcont,XAVE_sspcont_K,XSTD_sspcont,XSTD_sspcont_K,...
    YAVE_sspcont,YAVE_sspcont_T,YSTD_sspcont,YSTD_sspcont_T,XYCOR_sspcont,fignum};
end
%...................................................................................
%===================================================================================
%...................................................................................
if strcmp(key2T,'yes') && strcmp(keyDisc,'yes')
    fignum = [24];
    dat24 = {PHYTsspdisc,PHYTsspcont,logESDphysspAveDisc,logESDphysspAveCont,...
    TOPTphysspAveDisc,TOPTphysspAveCont,logESDphysspStdDisc,logESDphysspStdCont,...
    TOPTphysspStdDisc,TOPTphysspStdCont,physspCorDisc,physspCorCont,ndepths,fignum};
end
%...................................................................................

%===================================================================================
%...................................................................................
if strcmp(key2T,'yes') && strcmp(keyDisc,'yes')
    disp('Mean error and root mean square error on total phytoplankton, mean size and standard deviation (truth is the discrete model)')
    BiasPhy = mean(mean(PHYTsspcont - PHYTsspdisc)) / mean(mean(PHYTsspdisc))
    LSQPhy  = sqrt(mean(mean((PHYTsspcont - PHYTsspdisc).^2))) / mean(mean(PHYTsspdisc))
    BiasAve = mean(mean(logESDphysspAveCont - logESDphysspAveDisc)) % / mean(mean(logESDphysspAveDisc))
    LSQAve  = sqrt(mean(mean((logESDphysspAveCont - logESDphysspAveDisc).^2))) % / mean(mean(logESDphysspAveDisc))
    BiasStd = mean(mean(logESDphysspStdCont - logESDphysspStdDisc)) / mean(mean(logESDphysspStdDisc))
    LSQStd  = sqrt(mean(mean((logESDphysspStdCont - logESDphysspStdDisc).^2))) / mean(mean(logESDphysspStdDisc))
end
%...................................................................................
if strcmp(key2T,'yes')
    disp('Total primary production, mean phytoplankton, mean traits, variances and correlation weighted by phytoplankton in the continuous model')
    disp('Mnimum, average, and maximum')
    logDIN_minavemax = [min(log(DINsspcont(:))),sum(log(DINsspcont(:)).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(log(DINsspcont(:)))]
    T_minavemax = [min(temp(:)),sum(sum(temp(:,1:360).*PHYTsspcont))/sum(PHYTsspcont(:)),max(temp(:))]
    PP = deltaz*sum(PHYTsspcont(:).*MUPsspcont(:))*(106/16)*12*0.001 % in gC/m2/yr
    Phy_minavemax = [min(PHYTsspcont(:)),mean(PHYTsspcont(:)),max(PHYTsspcont(:))]
    logKnave_minavemax = [min(logESDphysspAveCont(:)),sum(logESDphysspAveCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(logESDphysspAveCont(:))]
    TOPTave_minavemax = [min(TOPTphysspAveCont(:)),sum(TOPTphysspAveCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(TOPTphysspAveCont(:))]
    logKnstd_minavemax = [min(logESDphysspStdCont(:)),sum(logESDphysspStdCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(logESDphysspStdCont(:))]
    TOPTstd_minavemax = [min(TOPTphysspStdCont(:)),sum(TOPTphysspStdCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(TOPTphysspStdCont(:))]
    Cor_minavemax = [min(physspCorCont(:)),sum(physspCorCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(physspCorCont(:))]
end
%...................................................................................
% Function to check convergence
if strcmp(key2T,'yes')
    [conv_P,conv_xave,conv_yave,conv_xxvar,conv_yyvar,conv_xycor] = SPEAD_1D_convergence_check(PHYTodecont,XAVE_odecont,YAVE_odecont,XXVAR_odecont,YYVAR_odecont,XYCOV_odecont,ndays);
end
%...................................................................................
%***********************************************************************************
% return

end %function []=SPEAD_1D_main()

