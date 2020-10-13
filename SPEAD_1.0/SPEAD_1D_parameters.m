function [deltat,nyear,ndays,deltaz,ndepths,betaz,gzmax,kgz,mz,mpower,Q10a,Q10h,temp0,Isat,InhFac,mup0,alp0,amup,aalp,aknp,mp,epsZoo,omePhy,omeZoo,md,galfa0,gbeta0,kw,wsink,xmin,xmax,xsigmaPcnt,ymin,ymax,ysigmaPcnt,nxphy,nyphy,numutX0,numutY0,phy0,zoo0,din0,pon0,box0] = SPEAD_1D_parameters(keyKTW,keyTraitDiffusion,keyModelResol)

%========================================================================
% Define all numerical, physical and biological constants of the model
%========================================================================

%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL CONSTANTS
%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................

% Time step (d)
% deltat = 1/48; %[days]
% deltat = 1/24; %[days] (ie. deltat = 1h)
% deltat = 2/24;
% deltat = 3/24; %[days] 
deltat = 6/24; %[days]
% deltat = 8/24; %[days]
% deltat = 12/24;
% deltat = 1; % Only possible if implicit scheme or step-splitting for vertical diffusion
%........................................................................
ndays = 360;
%........................................................................
%FOR REGULAR RUNS -- USE THREE YEARS OR MORE:
%%nyear = 1; 
%%nyear = 2; 
nyear = 3; %Standard
%%nyear = 10;
%........................................................................
%FOR FAST TESTING ONLY -- ** COMMENT OUT FOR REGULAR RUNS **
%%nyear = 1; 
%%ndays = ndays/12; %(ie. ndays = 30); %For fast-testing 0D runs.
%........................................................................
% Vertical step (m)
%%deltaz = 2.5; %[m] (needs dt < 0.015)
deltaz = 10.0; %[m] (needs dt < 0.25)
%%deltaz = 50.0; %[m]
%%deltaz = 100.0; %[m]
%........................................................................
if strcmp(keyModelResol,'0D')
    ndepths = 1; % 0D mode
elseif strcmp(keyModelResol,'1D')
    ndepths = 20; % 1D mode with 20 vertical levels
end
%........................................................................

%%%%%%%%%%%%%%%%%%%%%%
% BIOLOGICAL CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
mp = 0.05; % Phytoplankton mortality rate [d-1]
betaz = 0.4; % Zooplankton assimilation efficiency [%]
gzmax = 1.5; % Zooplankton maximum growth rate [d-1]
kgz = 0.4; % Zooplankton half-saturation [mmolN/m3]
mz = 0.25; % Zooplankton mortality rate [d-1]
mpower = 2.0; % Zoo mortality exponent (1=linear, 2=quadratic).
% &&& mpower = 1.0; %Linear mortality.
% &&& mpower = 1.4; 
%........................................................................
% Article values (24/08/2020)
Q10a = 1.75; % Eppley Q10 coefficient for autotrophic processes (1.0 = no temperature effect)
Q10h = 2.5; % Q10 coefficient for heterotrophic processes
temp0 = 20; % Reference temperature
%........................................................................
Isat = 25;
InhFac = 12;
%........................................................................
mup0 = 1.1; %Phy maximum growth rate [d-1] at Kn = 1 mmol/m3 and Topt = 20 deg
alp0 = mup0; % Affinity ((m3/mmol)/d) for Kn = 1 mmol/m3
amup = 0.5; % "Allometric" scaling of maximum growth rate (as a function of Kn)
aknp = 1; % By definition, since half-saturation is the trait
aalp = amup - aknp; % "Allometric" scaling of affinity
%........................................................................
% Standard values
%epsPhy = 1/3; %Phyplankton exudation fraction going to nitrate [n.d] 
epsZoo = 1/3; %Zooplankton exudation fraction going to nitrate [n.d]
omePhy = 1/4; %Phyplankton mortality fraction going to nitrate [n.d]
omeZoo = 1/4; %Zooplankton mortality fraction going to nitrate [n.d]
%........................................................................
md = 0.1; % Detritus remineralization rate at 20 degress [d-1]
%...................................................................................
if strcmp(keyKTW,'not')
    %...............................................................................
    galfa0 = 1.0; %no KTW.
    gbeta0 = 2.0; 
    %...............................................................................
elseif strcmp(keyKTW,'yes')
    %...............................................................................
% $$$     galfa0 = 1.2; %weak KTW.
% $$$     gbeta0 = 2.0; 
    %...............................................................................
    galfa0 = 1.5; %moderate KTW.
    gbeta0 = 2.0; 
    %...............................................................................
% $$$     galfa0 = 1.6; %moderate KTW.
% $$$     gbeta0 = 2.0; 
    %...............................................................................%
% $$$     galfa0 = 2.0; %strong KTW.
% $$$     gbeta0 = 2.0; 
    %...............................................................................
end
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%===================================================================================
%........................................................................
kw = 0.04; %irradiance attenutation due to water [m-1]
wsink = 1.2; % Detritus sinking speed [m/d]

%%%%%%%%%%%%%%%%%
% TRAIT CONSTANTS
%%%%%%%%%%%%%%%%%
%........................................................................
% Natural logarithm of half-saturation (in mmolN/m3)
xmin = -2.5;
xmax = 1.5;
xsigmaPcnt = 0.025;
%........................................................................
% Optimal temperature (degree)
ymin = 18;
ymax = 30;
ysigmaPcnt = 0.025;
%........................................................................
% Number of species in the discrete model
nxphy = 5; % Reference
nyphy = 5; % Reference
% Mutation rate per generation (for each trait)
if strcmp(keyTraitDiffusion,'yes')
    numutX0 = 0.003;
    numutY0 = 0.03;
elseif strcmp(keyTraitDiffusion,'not')
    numutX0 = 0;
    numutY0 = 0;
end

%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%
%........................................................................
phy0 = 0.1;
zoo0 = 0.1;
din0 = 1.8;
pon0 = 0.0;
box0 = 0.0;
%........................................................................

%........................................................................
%========================================================================
%************************************************************************
return


