function [deltat,nyear,ndays,deltaz,ndepths,parmin,parmax,skzmin,skzmax,kzdeep,sstmin,sstmax,betap_max,betap_xmax,betap_xrange,betaz,gzmax,kgz,mz,mpower,amuz,Q10a,Q10h,temp0,Isat,InhFac,mup0,alp0,amup,bmup,aalp,aknp,mp,epsPhy,epsZoo,omePhy,omeZoo,md,galfa0,gbeta0,kw,kp,wsink,ygamma,teta,xmin,xmax,xsigmaPcnt,ymin,ymax,ysigmaPcnt,nxphy,nyphy,numutX0,numutY0,phy0,zoo0,din0,pon0,box0] = jamstecrest_gaussecomodel1D_parameters(keyNutrientSupplyFlux,keyKTW,keyTraitDiffusion,keyModelResol)

%========================================================================
% Define all numerical, physical and biological constants of the model
%========================================================================

%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL CONSTANTS
%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
%Values by Sergio (24/08/2020)
%1 300 m/d
%20 120 w/m2 Cropp 120
%Q10 = 1
%nux = 0.05
%mphyto = 0.1 (background)
%mzoo =  0.2 (quadratic)

%mup0 = 1 (1 micron)
%alpha0 = 10
%amup = 0.5
%aalp = -amup
%bmup = 0;

%betaz = 0.4
%gzmax = 2.5
%kgz = 0.75
%amuz = 0.5 (?? -> 0)

% Time step (d)
% deltat = 1/48; %[days]
% deltat = 1/24; %[days] (ie. deltat = 1h) %USAR ESTE!!!
% deltat = 2/24;
% deltat = 3/24; %[days] 
deltat = 6/24; %[days]
% deltat = 8/24; %[days]
% deltat = 12/24;
% deltat = 1; % Only possible if implicit scheme or step-splitting for vertical diffusion (Le Gland, 23/10/2019)
%........................................................................
% Number of days in one year of simulation
if strcmp(keyNutrientSupplyFlux,'not')
    ndays = 360; %For months of 30 days (when using seasonal runs). 
elseif strcmp(keyNutrientSupplyFlux,'yes')
    ndays = 4*8*12; %For weeks of 8 days per month (when using double-period). 
end
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

%%%%%%%%%%%%%%%%%%%%
% PHYSICAL CONSTANTS
%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
parmin =  20.00; %winter surface PAR [W*m-2] %ORIGINAL VALUES.
parmax = 120.00; %summer surface PAR [W*m-2]
%........................................................................
% mldmin =  20.0; %min MLD [m]
% mldmax = 200.0; %max MLD [m]
%........................................................................
% Sergio's values (24/08/2020)
skzmin = 30;
skzmax = 300;
kzdeep = 1;
% 
%skzmin =   100; %min surface turbulence [m^2*day-1] %for "m4cMIT-v055.m" setup
%skzmax =  1000; %max surface turbulence [m^2*day-1]
%kzdeep =    10;
%........................................................................
%sstmin = 10; %min SST [degree]
%sstmax = 20; %max SST [degree]
sstmin = 12;
sstmax = 32;
%........................................................................

%%%%%%%%%%%%%%%%%%%%%%
% BIOLOGICAL CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
%betap_max = 0.8; %Phy assimilation efficiency [%]
betap_max = 1.0; % Maximum phytoplankton assimilation efficiency (Le Gland, 03/10/2019)
% betap_xmax and betap_xrange are used only if betap variability is activated
betap_xmax = 2.0; % Where assimilation efficiency is maximum (Le Gland, 03/10/2019)
betap_xrange = 3.0; % Width of Gaussian betap function (Le Gland, 03/10/2019)
% Phytoplankton mortality was not originally quadratic (Le Gland, 31/10/2019)
% Perhaps it should depend on temperature (heterotrophic or autotrophic Q10 ?)
%mp = 0.20; %Phy background mortality [d-1] %JAMSTEC original. (Quadratic)
% mp = 0.05; % Vallina et al. [2017] (Linear)
mp = 0.05; % Article value (24/08/2020)
% mp = 0.1; % Sergio's value (24/08/2020)
betaz = 0.4; %Zoo assimilation efficiency [%]
%gzmax = 1.0; %Zoo maximum grazing rate [d-1]
gzmax = 1.5; % Article value (24/08/2020)
%gzmax = 2.5; % Sergio's value (24/08/2020)
%gzmax = 2.0;
%kgz = 0.75; %Zoo half-sat grazing [mmolN*m-3] % Sergio's value (24/08/2020)
% kgz = 1.0; % Hansen et al. [1997] (Le Gland, 04/11/2019)
% kgz = 0.50;
kgz = 0.4; % Article value
%mz = 0.20; %Zoo background mortality [m3*mmolN-1*d-1] %JAMSTEC original. % Sergio's value (24/08/2020)
mz = 0.25;% Article value (24/08/2020)
mpower = 2.0; % Zoo (and phy ?) mortality exponent (1=linear, 2=quadratic).
% &&& mpower = 1.0; %Linear mortality.
% &&& mpower = 1.4; 
%amuz = -0.5; % Allometric grazing (and zooplankton growth rate)
amuz = 0;
%........................................................................
% $$$ betap = 1.0; %Phy assimilation efficiency [%]
% $$$ betaz = 0.8; %Zoo assimilation efficiency [%]
% $$$ gzmax = 1.0; %Zoo maximum grazing rate [d-1]
% $$$ kgz = 0.5; %Zoo half-sat grazing [mmolN*m-3]
% $$$ mz = 0.20; %Zoo background mortality [m3*mmolN-1*d-1] 
%........................................................................
% $$$ betap = 0.8; %Phy assimilation efficiency [%]
% $$$ betaz = 0.4; %Zoo assimilation efficiency [%]
% $$$ gzmax = 2.0; %Zoo maximum grazing rate [d-1]
% $$$ kgz = 0.50; %Zoo half-sat grazing [mmolN*m-3]
% $$$ mz = 0.20; %Zoo background mortality [m3*mmolN-1*d-1]
%........................................................................
% Sergio's values (24/08/2020)
%Q10a = 1.0;
%Q10h = 1.0;
% Article values (24/08/2020)
Q10a = 1.75; % Eppley Q10 coefficient for autotrophic processes (1.0 = no temperature effect)
Q10h = 2.5; % Q10 coefficient for heterotrophic processes (Le Gland, 31/10/2019)
temp0 = 20; % Reference temperature (Le Gland, 31/10/2019)
%........................................................................
% Isat = 120; %Phytoplankton saturating irradiance [W*m-2] % Sergio's value (24/08/2020)
Isat = 25; % I think this is closer to Follows et al. [2007] and it allows for growth of deep plankton (Le Gland, 31/10/2019) % Article value (24/08/2020)
%Isat = 5/(log(sqrt(3))); % correspond to a half-saturation of 5 in a hyperbolic tangent model
% InhFac = 12; % Photoinhibition factor (KPAR/Kinhib in Follows et al. [2007])
% InhFac = 0.5 makes it very similar to the formulation by Cropp (2004) and Walsh (2001) used previously
InhFac = 12;
%........................................................................
%mup0 = 1.00; %Phy maximum grazing rate [d-1] at ESD = 1.0 %ORIGINAL VALUES.
%mup0 = 1.50; % To account for non-linear maximum growth rate (Le Gland, 25/09/2019)
%alp0 = 10.0; %Phy maximum uptake affinity at ESD = 1.0 [d-1.mmolN-1.m3]
%alp0 = 15.0; % To account for non-linear maximum growth rate (Le Gland, 25/09/2019)
%mup0 = 0.5; % More realistic, avoids cells of 30 microns (the largest possible here) to grow larger than the Eppley curve (Le Gland, 25/10/2019)
%alp0 = 5.0; % Species of 1 micron must have a Kd for N of 0.1 mmolN.m-3 and Kd grows like size (Edwards et al., 2012) (Le Gland, 25/10/2019)
% Case where knp is the x trait
%mup0 = 1.5; %Phy maximum growth rate [d-1] at Kn = 1 mmol/m3 and T = 20 deg
%mup0 = 1.2;
mup0 = 1.1; % Article value (24/08/2020)
%mup0 = 3.1; % Sergio's value (24/08/2020)
alp0 = mup0; % By definition
%amup = 1.0; %Size scaling factor for muPhy [n.d.] JAMSTEC original.
amup = 0.5; % More realistic value ? (Le Gland, 20/09/2019)
%amup = 0.389; % Value for variable assimilation efficiency (Le Gland, 03/10/2019)
bmup = 0;
%bmup = -0.055; %Quadratic scaling factor for muPhy
% aalp = -amup; %Size scaling factor for alPhy [n.d.]
% aknp = (amup - aalp); %Size scaling factor for ksPhy [n.d.]
% Adapted from Chen and Smith (2018)
%amup = 0.6;
%bmup = -0.03;
%aalp = amup - 0.8;
%aknp = (amup - aalp); % 0.8 (Litchman et al., 2007)
% Case where knp is the x trait (Le Gland, 31/10/2019)
aknp = 1; % By definition
aalp = amup - aknp;
%........................................................................
% Standard values
epsPhy = 1/3; %Phyplankton exudation fraction going to nitrate [n.d] 
epsZoo = 1/3; %Zooplankton exudation fraction going to nitrate [n.d]
omePhy = 1/4; %Phyplankton mortality fraction going to nitrate [n.d]
omeZoo = 1/4; %Zooplankton mortality fraction going to nitrate [n.d]
%........................................................................
% md = 0.10; %Detritus degradation specific rate [d-1] (original).
% Should increase with temperature with the Q10 of heterotrophs, but what
% rate at 20 deg ? (Le Gland, 31/10/2019)
% md = 0.15; % 0.1 at 15 deg, reference temp in Chen et al., 2018
% md = 0.08;
md = 0.1;
% $$$ md = 2.00; %Detritus degradation specific rate [d-1] (super fast recycling)
%........................................................................
%TO REMOVE DETRITUS FROM THE SYSTEM MAKING IT A NPZ MODEL:
%........................................................................
% epsPhy = 1; %Phyplankton exudation fraction going to nitrate [n.d] 
% epsZoo = 1; %Zooplankton exudation fraction going to nitrate [n.d]
%..................................................................
% omePhy = 1; %Phyplankton mortality fraction going to nitrate [n.d]
% omeZoo = 1; %Zooplankton mortality fraction going to nitrate [n.d]
%..................................................................
% md = 0.00; %Detritus degradation specific rate [d-1]
%........................................................................
%TO REMOVE NUTRIENT RECYCLING FROM THE SYSTEM MAKING: 
%........................................................................
% $$$ epsPhy = 0; %Phyplankton exudation fraction going to nitrate [n.d] 
% $$$ epsZoo = 0; %Zooplankton exudation fraction going to nitrate [n.d] 
% $$$ %..................................................................
% $$$ omePhy = 0; %Phyplankton mortality fraction going to nitrate [n.d] 
% $$$ omeZoo = 0; %Zooplankton mortality fraction going to nitrate [n.d] 
% $$$ %..................................................................
% $$$ md = 0.00; %Detritus degradation specific rate [d-1] 
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
kp = 0.00; %irradiance attenuation due to phyto self-sheding [m2*mmolN-1]
% kp = 0.02; % Rate from Vallina et al. (2017) (Le Gland, 04/11/2019)
% Should depend on Chl:N ratio, which is variable (Le Gland, 04/07/2019)
% wsink = 1.0; % Particle sinking speed [m*d-1]
% wsink = 1.6;
wsink = 1.2; % Article value (24/08/2020)
%........................................................................
ygamma = 4.0; % Width of temperature window for a given species
teta = 2;

%%%%%%%%%%%%%%%%%
% TRAIT CONSTANTS
%%%%%%%%%%%%%%%%%
%........................................................................
% Natural logarithm of Equivalent Spherical Diameter (ESD, in um)
% xmin = -3;
% xmax = 5;
% Replace ESD by Half-saturation (mmol/m3) (Le Gland, 31/10/2019)
% xmin = -2.9957; % 0.05 mmol/m2
% xmax = 1.6094; % 5 mmol/m3
xmin = -2.5;
xmax = 1.5;
% xmin = -1.0;
% xmax = 1.0;
%xsigmaPcnt = 0.1;
%xsigmaPcnt = 0.05;
xsigmaPcnt = 0.025;
%........................................................................
% Optimal temperature (degree)
ymin = 18;
ymax = 30;
% ymin = 21;
% ymax = 27;
%ysigmaPcnt = 0.1;
%ysigmaPcnt = 0.05;
ysigmaPcnt = 0.025;
%........................................................................
% Optimal irradiance (when available as a trait)
% Number of species in the discrete model
nxphy = 25; % Reference
nyphy = 25; % Reference
%nxphy = 3; % Number of values for the x trait
%nyphy = 3; % Number of values for the y trait
% Trait diffusivity parameter (for each trait)
if strcmp(keyTraitDiffusion,'yes')
    % Reach the maximum variance ( ((max-min)/2)^2 ) in 1000 generations
    numutX0 = 0.001;
    numutY0 = 0.01;
elseif strcmp(keyTraitDiffusion,'not')
    numutX0 = 0;
    numutY0 = 0;
end

%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%
%........................................................................
%OLIGOTROPHIC:
%phy0 = 1.0; %Phy [mmolN*m-3]
%zoo0 = 1.0; %Zoo [mmolN*m-3]
%din0 = 1.0; %DIN [mmolN*m-3]
%pon0 = 0.0; %PON [mmolN*m-3]
%box0 = 0.0; %BOX [mmolN*m-3]
%...................................................................................
%OLIGOTROPHIC 2
% $$$ phy0 = 0.1; %Phy [mmolN*m-3]
% $$$ zoo0 = 0.1; %Zoo [mmolN*m-3]
% $$$ din0 = 2.8; %DIN [mmolN*m-3]
% $$$ pon0 = 0.0; %PON [mmolN*m-3]
% $$$ box0 = 0.0; %BOX [mmolN*m-3]
%........................................................................
%EUTROPHIC:
% $$$ phy0 = 1.0; %Phy [mmolN*m-3]
% $$$ zoo0 = 1.0; %Zoo [mmolN*m-3]
% $$$ din0 = 4.0; %DIN [mmolN*m-3]
% $$$ pon0 = 0.0; %PON [mmolN*m-3]
% $$$ box0 = 0.0; %BOX [mmolN*m-3]
%MORE OLIGOTROPHIC (Le Gland, 16/12/2019)
phy0 = 0.1;
zoo0 = 0.1;
din0 = 1.8; % Article value (24/08/2020)
%din0 = 2.8; % Sergio's value (24/08/2020)
pon0 = 0.0;
box0 = 0.0;
%........................................................................

%........................................................................
%========================================================================
%************************************************************************
return


