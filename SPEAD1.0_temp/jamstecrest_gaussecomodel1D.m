%********************************************************************
% PROGRAM: JAMSTEC_GAUSSCOMODEL1D.M
%********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FILE HEADER -- LOAD PACKAGES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************************************************************

clear all
%--------------------------------------------------------------------
addpath('../CODE-RUNS/MYFUNCTIONS/');
%--------------------------------------------------------------------
addpath(genpath('../CODE-RUNS/MYTOOLBOX/MATLAB-LIBRERY/'));
%--------------------------------------------------------------------
%====================================================================
%--------------------------------------------------------------------
[mypackages] = myheadloadpackages; %Structure array to pass on my head pkg as input argument to functions.
%--------------------------------------------------------------------
%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal;
%--------------------------------------------------------------------
%JUST CHECKING IF PLOTING WORKS OKAY:
A256 = peaks(256);
A128x256 = A256(1:2:256,:);
A = A128x256;
Amin = min(A(:));
Amax = max(A(:));
fignum = 1001;
[hcbar] = jamstecrest_subplotesting(A,Amin,Amax,fignum,mypackages);
pause(0.5)
close all
%--------------------------------------------------------------------
%====================================================================
%********************************************************************
global galfa gbeta
% global gzmax kgz mz betaz betap mpower 
global gzmax kgz mz betaz betap_max betap_xmax betap_xrange mpower 
global mp Isat InhFac numutx numuty % Replace numut by 2 different mutation rates (Le Gland, 10/09/2019) 
global alp0 mup0 % knp0 
global aalp amup aknp bmup 
global amuz % Zooplankton dependent grazing (Le Gland, 26/08/2020)
%global Q10 
global Q10a Q10h % Distinct partition coefficients for auto and heterotrophic processes (Le Gland, 31/10/2019)
global temp0 sst temp % sst0 Use temperature at all depths (Le Gland, 15/10/2019)
global tcounter 
global jcounter 
%global icounter 
%global iTimeode 
global keyTraitAxis keyPhysics keySinking keyAssimConstant % keyAssimConstant added by Le Gland, 03/10/2019
global Sdin Mdin Drate 
global epsPhy omePhy epsZoo omeZoo md 
global t0 deltat ndays nyear tmax tspan 
global zdepths ndepths deltaz
global Ixave Iyave Ixxvar Iyyvar Ixycov % continuous model with 2 traits (Le Gland, 02/07/2019)
global Ixave_K Ixxvar_K Iyave_T Iyyvar_T % 1-trait model indices (Le Gland, 25/11/2019)
global Iphy Izoo Idin Ipon Ibox %continuous model
global Jphy Jzoo Jdin Jpon Jbox %discrete model  
global nxphy nyphy nzoo ndin npon nbox % discrete model with 2 traits (Le Gland, 04/07/2019)
global ntot0 %continuous model
global mtot0 %discrete model
global KZ % KZI
global parz0
global kw kp wsink
global jday jjday 
global xtrait ytrait 
global myTitle001 myTitle002 myTitle003 myTitle004 
global myTitle221 myTitle222 myTitle223 myTitle224 myTitle225 myTitle226 % myTitle225 added by Le Gland (05/09/2019)
global myXtickMarks myXtickLabel myXaxisLabel
global myYtickMarks myYtickLabel myYaxisLabel
%global keyFastNumericalSolving
global keyNutrientSupplyFlux 
global keyTraitDiffusion
global xrng yrng dx dy % 2 traits (Le Gland, 16/07/2019)
global KPPXstar KPPXstarI 
%...................................................................................
global Xavedotday Yavedotday XXvardotday YYvardotday XYcovdotday %OUTPUTS (2 traits) 
global UXYday GXYday 
global d1UXYdxday d1UXYdyday d1GXYdxday d1GXYdyday 
global d2UXYdxdxday d2UXYdydyday d2UXYdxdyday d2GXYdxdxday d2GXYdydyday d2GXYdxdyday
global todedotday
global Xavedotout Yavedotout XXvardotout YYvardotout XYcovdotout %OUTPUTS (2 traits)
global UXYout GXYout
global UXout GXout UYout GYout % 1-trait model outputs (Le Gland, 27/11/2019)
global d1UXYdxout d1UXYdyout d1GXYdxout d1GXYdyout 
global d2UXYdxdxout d2UXYdydyout d2UXYdxdyout d2GXYdxdxout d2GXYdydyout d2GXYdxdyout
global todedotout
%...................................................................................
global uphydaydisc gphydaydisc %OUTPUTS FROM ODE45 OF DISCRETE MODEL. 
global uphyoutdisc gphyoutdisc %OUTPUTS FROM ODE45 OF DISCRETE MODEL. 
%...................................................................................
global FPHYdaydisc GPHYdaydisc %OUTPUTS FROM ODE45 OF DISCRETE MODEL. 
global FPHYoutdisc GPHYoutdisc %OUTPUTS FROM ODE45 OF DISCRETE MODEL. 
%...................................................................................
global FPHYToutcont EPHYToutcont MPHYToutcont GPHYToutcont %OUTPUTS FROM ODE45 OF CONTINOUS MODEL. 
global FZOOoutcont  EZOOoutcont  MZOOoutcont 
global FDINoutcont  FPONoutcont
%...................................................................................
global FPHYToutdisc EPHYToutdisc MPHYToutdisc GPHYToutdisc %OUTPUTS FROM ODE45 OF DISCRETE MODEL. 
global FZOOoutdisc  EZOOoutdisc  MZOOoutdisc 
global FDINoutdisc  FPONoutdisc
%...................................................................................
global PHYmax ZOOmax DINmax PONmax 
global PHYmin ZOOmin DINmin PONmin 
%...................................................................................
global logESDaveMax logESDstdMax 
global logESDaveMin logESDstdMin 
global    ESDaveMax    ESDstdMax 
global    ESDaveMin    ESDstdMin 
global   TOPTaveMax   TOPTstdMax % Le Gland, 02/09/2019
global   TOPTaveMin   TOPTstdMin % Le Gland, 02/09/2019
global         CorrelationAbsMax % Le Gland, 05/11/2019 
global   EntropyMin   EntropyMax % Le Gland, 05/11/2019
%...................................................................................
global UXmin UXmax
global GXmin GXmax 
%...................................................................................
global MUPmin MUPmax
global MUZmin MUZmax 
%...................................................................................
global minArray2D maxArray2D 
%...................................................................................
global DIFFxaveout DIFFxvarout DIFFxstdout 
%...................................................................................
global myYlims1m myYtick1m myYtickLabel1m 
global myXlims1m myXtick1m myXtickLabel1m
global monthlimits
%...................................................................................

tic 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL KEYS AND PARAMETERS PACKAGES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%KEYS
%........................................................................
[keyKN,keyTOPT,key2T,keyModelResol,keyPhysics,keySinking,keyPARseasonality,...
keyPARextinction,keyNutrientSupplyFlux,keyNutrientPulses,keyKTW,keyTraitDiffusion,...
keyGalfaConstant,keyTraitDiffusionConstant,keyNutrientSupplyFrequencyConstant,keyAssimConstant] = ...
jamstecrest_gaussecomodel1D_keys;
%........................................................................
%========================================================================
%PARAMETERS
%........................................................................
[deltat,nyear,ndays,deltaz,ndepths,parmin,parmax,skzmin,skzmax,kzdeep,sstmin,sstmax,...
betap_max,betap_xmax,betap_xrange,betaz,gzmax,kgz,mz,mpower,amuz,Q10a,Q10h,temp0,Isat,InhFac,mup0,alp0,amup,bmup,aalp,aknp,mp,epsPhy,epsZoo,omePhy,omeZoo,...
md,galfa0,gbeta0,kw,kp,wsink,ygamma,teta,xmin,xmax,xsigmaPcnt,ymin,ymax,ysigmaPcnt,nxphy,nyphy,numutX0,numutY0,...
phy0,zoo0,din0,pon0,box0] = jamstecrest_gaussecomodel1D_parameters(keyNutrientSupplyFlux,keyKTW,keyTraitDiffusion,keyModelResol);
%........................................................................

%%%%%%%%%%%%%%%%%%%%%
%TEMPORAL RESOLUTION:
%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
t0 = deltat;
tmax = ndays*nyear; 
tspan = [t0:deltat:tmax]; 
%........................................................................
nsteps = length(tspan); % Should be equal to (ndays/deltat)*nyear
ttmax = ndays; %[days]
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%
%VERTICAL RESOLUTION:
%%%%%%%%%%%%%%%%%%%%%
%========================================================================

% Change in the depth levels (Le Gland, 06/09/2019)
% ndepths is now defined in the parameters
% First level begins at deltaz/2 (I view it as a box from 0 to deltaz)
zmin = 0;
zmax = ndepths*deltaz;
zdepths = [zmin+deltaz/2:deltaz:zmax-deltaz/2];

%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%
%EXTERNAL FORCINGS:
%%%%%%%%%%%%%%%%%%%

% Temperature is from a monthly climatology (arithmetic mean + interpolation with interp), Kz is from model GOTM (geometric mean)
% and PAR is from model GOTM (Le Gland, 14/10/2019)
% Should be re-interpolated if the depth step is changed
%load('../INPUTS/inputs.mat')
%iparz0 = parz0_gotm;
% iKZ = [ikzgotminterp(1,:);ikzgotminterp;ikzgotminterp(end,:)]; % In this file, KZ has been divided by 100 compared to GOTM estimates
% iKZ = [ikzgotminterp(1,:);ikzgotminterp;ikzgotminterp(end,:)]*100;% 
% isst = itempinterp(1,:);
% itemp = itempinterp;

% Load observations of PP, Chl, NO3 and PON (Le Gland, 05/12/2019)
% Units to be checked (Le Gland, 05/12/2019)
CHL_obs = load('../INPUTS/CHLbatsClimMonthlyProfilesPolyregress.txt');
PP_obs  = load('../INPUTS/PPbatsClimMonthlyProfilesPolyregress.txt');
PP_obs  = PP_obs * (1/12) * (16/106); % Converts mgC/m3/day to mmol/m2/day
NO3_obs = load('../INPUTS/NO3batsClimMonthlyProfilesPolyregress.txt');
PON_obs = load('../INPUTS/PONbatsClimMonthlyProfilesPolyregress.txt');
PON_obs = PON_obs / 14; % Converts ug/kg to mmol/m3

% GOTM is NOT in the Sargasso Sea. Now all forcings are from BATS raw measurements (Le Gland, 24/10/2019)
load('../INPUTS/MLDlevitus360days.txt')
load('../INPUTS/PARbats360days.txt')
load('../INPUTS/TEMPbatsinterp2D.txt')

% Interpolate cliamotological temperature into a ndepths*ndays forcing file
ndepths_temp = size(TEMPbatsinterp2D,1); % Number of depth level in temperature input file (200)
width_temp = floor(ndepths_temp / ndepths); % Ratio of depth levels in the input and forcing files
itemp = zeros(ndepths,ndays);
for i = 1:ndepths
    itemp_month = mean(TEMPbatsinterp2D(width_temp*(i-1)+1:width_temp*i,:),1);
    itemp(i,:) = interp1(1:365,itemp_month,1:364/359:365); % Interpolates on 360 days instead of 365
end
isst = itemp(1,:);

%NETCDF cannot be read from Octave (Le Gland, 15/07/2020)
% Read Diffusivity in a NETCDF file (from GOTM model)
% I could do it this way but I will not, since it is very noisy
% I will use their values as orders of magnitude (Le Gland, 05/12/2019)
% Perhaps I can use it, taking the average of the three years and
% smoothing with a time scale of 30 days (Le Gland, 11/12/2019)
ncid  = netcdf.open('../INPUTS/phys_BATS_for_Sergio.nc','NOWRITE');
varid = netcdf.inqVarID(ncid,'nuh');
Kz1   = netcdf.getVar(ncid,varid);
Kz2   = 60*60*24*squeeze(Kz1); % 1x1x250x1095 (m2/s) --> 250x1095 (m2/d)
Kz3   = flipud(Kz2(51:250,:));
Kz4   = (Kz3(:,1:365) + Kz3(:,366:730) + Kz3(:,731:1095)) / 3; % Average of three years (200x365) 
netcdf.close(ncid);

% Interpolate GOTM Kz into a ndepths*ndays forcing file and smooth
ndepths_Kz = size(Kz4,1); % Number of depth level in temperature input file (200)
width_Kz   = floor(ndepths_Kz / ndepths); % Ratio of depth levels in the input and forcing files
iKZ = zeros(ndepths+1,ndays);
%Kz5 = zeros(ndepths,365);
%Kz6 = zeros(ndepths,ndays);
%Kz7 = zeros(ndepths,390);
%Kz8 = zeros(ndepths,390);
for i = 1:ndepths+1
    Kz5 = mean(log(Kz4(max(1,width_Kz*(i-3/2)+1):min(width_Kz*(i-1/2),200),:)),1); % Interpolation in log-scale
    Kz6 = interp1(1:365,Kz5,1:364/359:365); % Interpolates on 360 days instead of 365
    Kz7 = [Kz6(346:360),Kz6,Kz6(1:15)];
    Kz8 = smooth(Kz7,30,'sgolay');
    iKZ(i,:) = exp(Kz8(16:375));
end

% Approximate transformation of Einstein/m2/d into W/m2
iparz0 = 2.5*PARbats360days; % Actually, depends on wavelength

% Make sure that the iparz0 max is at day 171, and not at day 180 as in
% PARbats360days (Le Gland, 15/01/2020)
iparz0 = zeros(1,360);
iparz0(1:351)   = 2.5*PARbats360days(10:360);
iparz0(352:360) = 2.5*PARbats360days(1:9); 

% Moved here to be plotted in the forcing figure (Le Gland, 31/10/2019)
PAR2D = exp(-kw*zdepths(:))*iparz0; % PAR profiles [W*m-2] (depth,time)

% Deduce diffusion coefficients (KZ) from mixed layer depth (MLD)
% 10000 m2/d at maximum MLD, 1000 at minimum, 10 below the MLD
imld = MLDlevitus360days;
mldmin = min(imld);
mldmax = max(imld);

%skzmax  = 1200;
%skzmin  = 60;
%kzdeep  = 3;
iskz = skzmax - (skzmax-skzmin)*((mldmax-imld)./(mldmax-mldmin));
%iKZ = jamstecrest_KZdiffprofiles(ndays,zmax,deltaz,imld,iskz,kzdeep);

fignum = 14;
jamstecrest_imagescforcings(itemp,iparz0,PAR2D,imld,iKZ,fignum,mypackages)

%=======================================
%PHOTOSYNTHETIC ACTIVE RADIATION (PAR+):
%=======================================
%........................................................................
% [iparz0] = jamstecrest_SinusoidalFunction(parmin,parmax,ndays,ttmax,'Linear'); %Seasonal PAR.
% [iparz0] = jamstecrest_SinusoidalFunction(parmin,parmax,ndays,ttmax,0,'Linear'); %Seasonal PAR.
%........................................................................

%========================
%MIXED LAYER DEPTH (MLD):
%========================
%........................................................................
% MLD and temperature are shifted by 2 months compared to PAR (Le Gland, 27/09/2019)
% Not sure if it is realistic for MLD
% shift = ndays/6;  
% [imld] = jamstecrest_SinusoidalFunction(mldmax,mldmin,ndays,ttmax,'Quadratic');
% [imld] = jamstecrest_SinusoidalFunction(mldmax,mldmin,ndays,ttmax,shift,'Quadratic');
%........................................................................

%===================================
%SURFACE TURBULENCE DIFFUSION (SKZ):
%===================================
%........................................................................
% iskz = skzmax - (skzmax-skzmin)*((mldmax-imld)./(mldmax-mldmin)); %skz is a function of the MLD.
%........................................................................

% I now try to use only one set of KZ, with values taken at the boundaries of boxes
% No interpolation will be required any longer (Le Gland, 06/09/2019)
% if ndepths > 1 % There is no vertical diffusion in a 0D model
%     % Turbulent diffusion at nodes j-1/2 and j+1/2
%     iKZ = jamstecrest_KZdiffprofiles(ndays,zmax,deltaz,imld,iskz);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SINUSOIDAL PULSED FUNCTION OF NUTRIENT SUPPLY (S):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
if strcmp(keyNutrientSupplyFlux,'not')
    %...............................................................................
    %HAVING ZERO DIN FLUX AND ZERO OVERFLOW DILUTION:
    %...............................................................................
    Drate0 =  0.0; %Constant dilution specific rate [d-1]
    DINext0 = 0.0; %Concentration of DIN exterior to my volume of water [mmolN*m-3]
    %...............................................................................
    sdin  = Drate0 * DINext0; %Flux of nutrient supply [mmolN*m-3*d-1]
    drate = Drate0; %[d-1] 
    %...............................................................................
    jSdin  = sdin  * ones(1,ndays/deltat); 
    jDrate = drate * ones(1,ndays/deltat); 
    %...............................................................................
elseif strcmp(keyNutrientSupplyFlux,'yes')
    %...............................................................................
    nfreqs = 8; 
    %%nfreqs = ndepths;
    %...............................................................................
    [SDIN,DRATE,fsin,ttime,Wrad] = jamstecrest_funsin1D(nfreqs,ndays,t0,deltat,keyNutrientPulses,keyNutrientSupplyFrequencyConstant);
    %...............................................................................
    if exist('jloop') == 0 %if does not exist. 
	    jloop = nfreqs;  
    end
    %...............................................................................
    jSdin = SDIN(jloop,:); 
    jDrate = DRATE(jloop,:); 
    %...............................................................................
end
%...................................................................................
[msize,nsize] = size(jSdin);
%...................................................................................
if msize == ndepths
    iSdin  = jSdin;
    iDrate = jDrate;
elseif msize == 1; 
    iSdin  = ones(ndepths,1) * jSdin; 
    iDrate = ones(ndepths,1) * jDrate; 
end
%...................................................................................
%===================================================================================
%%return

%%%%%%%%%%%%%%%%%%%%%%%%%
%SEA SURFACE TEMPERATURE:
%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
% [isst] = jamstecrest_SinusoidalFunction(sstmin,sstmax,ndays,ttmax,'Linear');
% [isst] = jamstecrest_SinusoidalFunction(sstmin,sstmax,ndays,ttmax,shift,'Linear');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAKE LONGER ARRAYS (mld, par, KZ, KZI) ACCORDING TO NUMBER OF YEARS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%........................................................................
if ndepths > 1
    KZ  = repmat(iKZ,[1,nyear]);
    % KZI = repmat(iKZI,[1,nyear]); % No longer necessary (Le Gland, 06/09/2019)
end
mld = repmat(imld,[1,nyear]);
sst = repmat(isst,[1,nyear]);
temp = repmat(itemp,[1,nyear]); % Use temperature at all depths (Le Gland, 15/10/2019)
Sdin = repmat(iSdin,[1,nyear]);
Drate = repmat(iDrate,[1,nyear]);
parz0 = repmat(iparz0,[1,nyear]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STABILITY CONDITION FOR TURBULENCE SCHEME: 
%(not really relevant when using "ode45" solver 
% because it uses a variable adaptive dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------
%NOTE: Algo esta MAL!!! con el sinking alone (check it out)
%----------------------------------------------------------
%======================================
%--------------------------------------
% dtdiff < min(dz^2/(2*kzI))
%--------------------------------------
%........................................................................
%wsink = 1.0; %[m*d-1]
%........................................................................
if ndepths > 1
    %maxKZ = max(KZI);
    maxKZ = max(KZ); % KZ is no longer interpolated (Le Gland, 06/09/2019)
    stability = min(deltaz^2./(2*maxKZ)) %This defines the maximum dt allowed.
    %........................................................................
    if (deltat*2) >= stability 
        disp('Be careful: Condition of stability violated!')
        disp('reduce dt or increase dz')
        pause
    end
    %........................................................................
    dtadv = min((deltaz^2)./((wsink*deltaz) + (2*maxKZ))); %max time step from von Neumann stability analysis
    %........................................................................
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHY CELL SIZE RESOLUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%if strcmp(keyTraitAxis,'ESD')
%...................................................................................
    % Define ESDmin and ESDmax based on the x trait (Le Gland, 09/09/2019)
    ESDmin = exp(xmin);
    ESDmax = exp(xmax); 
    xdel = (xmax-xmin)/(nxphy-1); % Probably the most consistent expression, since no point has to be removed (Le Gland, 04/06/2019)   
    %...............................................................................
    xaxis  = [xmin:xdel:xmax]; %[log(um)] 
    %...............................................................................
	ESD = exp(xaxis); %Equivalent Spherical Diameter [um]
    %...............................................................................
    show_ESD = [[1:nxphy]',xaxis(:),ESD(:)]
%elseif strcmp(keyTraitAxis,'SST')
    ydel = (ymax-ymin)/(nyphy-1);
    yaxis = [ymin:ydel:ymax];
%end
%...................................................................................
%===================================================================================

%...................................................................................
%===================================================================================
%TRAIT DIFFUSION DISPERSION: 
%...................................................................................
% Direct expression of trait diffusivity parameter (Le Gland, 01/07/2019)
% sigmaXm and sigmaYm are OK for analyses, but are not necessary in the 
% computations (Le Gland, 18/10/2019)
sigmaXm = sqrt(numutX0);
sigmaYm = sqrt(numutY0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IF YOU WANT TO REMOVE ALL PHYSICAL PROCESSES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
if strcmp(keyPhysics,'not') && ndepths > 1
    KZ (:,:) = 0; %no turbulence
    % KZI(:,:) = 0; %no turbulence % KZI is no longer required (Le Gland, 09/09/2019)
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
end
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTAIN 2D FIELDS FOR PLOTTING IMAGESC LATER:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
% PAR2D = exp(-kw*zdepths(:))*iparz0; % PAR profiles [W*m-2] (depth,time) (Le Gland, 01/07/2019)
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%
%GAUSSIAN DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%
%if strcmp(keyTraitAxis,'ESD')
    %------------------------------------------------------------------------
    %NOTE: fxj = (1.0 ./ (sigmaxj * sqrt(2*pi))) .* exp( -(xj - xmj).^2 ./ (2*sigmaxj.^2) ); %Okay.
    %------------------------------------------------------------------------
    %===================================================================================
    %-----------------------------------------------------------------------------------
    %NOTE: mean(x) = log(median(ESD))
    %-----------------------------------------------------------------------------------
    %...................................................................................
    %ESDave = exp(mean([xmin:xmax])); 
    ESDave = exp(mean([xmin,xmax])); % This is the correct mean (Le Gland, 01/07/2019) 
    %%ESDsig = (1/3)*ESDave; 
    %...................................................................................
    % $$$ % $$$ xmean = log(ESDave); 
    % $$$ % $$$ %%xsigma = log(ESDsig); %WRONG!!! (can lead to negative values!!!)
    %...................................................................................
    xmean = log(ESDave);
    %xsigmaPcnt = (1/10);  % normalized standard deviation (Le Gland, 29/08/2019)
    xsigma = xsigmaPcnt*(xmax-xmin); %Okay. 
    %...................................................................................
    fx = (1.0 / (xsigma * sqrt(2*pi))) * exp( -(xaxis - xmean).^2 / (2*xsigma^2) ); %Okay.
    %...................................................................................
    sumfx = sum(fx * xdel); %must be equal to one!!! 
    %...................................................................................
    if abs(sumfx - 1.0) > 1d-6 %test.
        disp('Error!!! sumfx must be equal to one!!!')
    end
    %...................................................................................
    xplus = xmean + xsigma;
    xless = xmean - xsigma;
    %...................................................................................
    xdelta = xplus - xless; 
    %...................................................................................
    ESDplus = exp(xmean) * exp(xsigma);
    ESDless = exp(xmean) / exp(xsigma);
    %...................................................................................
    ESDmin = exp(xmin);
    ESDmax = exp(xmax);
    %...................................................................................
    %===================================================================================
    % $$$ %...................................................................................
    % $$$ ESDvec = exp(xaxis);
    % $$$ ESDave = exp(xmean); %median!!!! (*not* mean!)
    % $$$ ESDsigma = exp(xsigma); %deviation 
    % $$$ %...................................................................................
    % $$$ fESD = (1.0 / (ESDsigma * sqrt(2*pi))) * exp( -(ESDvec - ESDave).^2 / (2 * ESDsigma^2) ); %ESTA TODO MAL!!!!
    % $$$ %...................................................................................
    % $$$ ESDmsigplus = ESDave * ESDsigma;
    % $$$ ESDmsigless = ESDave / ESDsigma;
    % $$$ %...................................................................................
    % $$$ ESDmsigdelta = ESDmsigplus - ESDmsigless; 
    % $$$ ESDmsigdeltaBis = (exp(xmean) / exp(xsigma)) * (exp(xsigma)^2 - 1); 
    % $$$ %...................................................................................
    %===================================================================================
    %-----------------------------------------------------------------------------------
    % <http://www.math.uah.edu/stat/special/Normal.html>
    % <http://www.math.uah.edu/stat/special/LogNormal.html>
    %-----------------------------------------------------------------------------------
    %...................................................................................
    % M = log(ESDave); %median of the distribution in ESD linear scale (mean in log(ESD) scale).
    % S = xsigmaPcnt*(log(ESDmax) - log(ESDmin));
    %...................................................................................
    % flogESD001 = (1.0  / (S * sqrt(2*pi)  ))    * exp( - (log(ESD) - M).^2 / (2 * S^2) ); %Wrong... (but actually nicer; I dont understand)
    % flogESD002 = (1.0 ./ (S * sqrt(2*pi)*ESD)) .* exp( - (log(ESD) - M).^2 / (2 * S^2) ); %Okay...?
    % Trait (logsize) is normally distributed, so flogESD001 is the right entropy (Le Gland, 01/07/2019)
    %...................................................................................
    % flogESD = flogESD001; %Seems okay!!!!
    % $$$ flogESD = flogESD002; %Seems *not* okay (something is going on; I dont understand)
    % Use xmean and xsigma instead of defining new variables (Le Gland, 29/08/2019)
    flogESD = (1.0  / (xsigma * sqrt(2*pi)  ))    * exp( - (xaxis - xmean).^2 / (2 * xsigma^2) );
    %...................................................................................
    % ESDplus = exp(M) * exp(S);
    % ESDless = exp(M) / exp(S);
    %...................................................................................
    ESDdelta = ESDplus - ESDless;
    %...................................................................................
    %===================================================================================
    figure(1000)
    %...................................................................................
    subplot(2,2,1)
    plot(xaxis,fx,'k-',xaxis,fx,'b.')
    hold on
    plot(xless,0,'r*',xplus,0,'r*')
    hold off
    set(gca,'Xlim',[xmin xmax])
    set(gca,'Ylim',[0.00 1.00])
    xlabel('log (size)')
    ylabel('f (x)')
    grid on
    %...................................................................................
    % $$$ subplot(2,2,2)
    % $$$ plot(ESDvec,fESD,'k-',ESDvec,fESD,'b.')
    % $$$ hold on
    % $$$ plot(ESDmsigless,0,'r*',ESDmsigplus,0,'r*')
    % $$$ hold off
    % $$$ set(gca,'Xlim',[exp(xmin) exp(xmax)])
    % $$$ set(gca,'Ylim',[0.00 1.00])
    % $$$ xlabel('(size)')
    % $$$ ylabel('f (ESD)')
    % $$$ grid on
    %...................................................................................
    subplot(2,2,3)
    plot(ESD,fx,'k-',ESD,fx,'b.')
    hold on
    plot(ESDless,0,'r*',ESDplus,0,'r*')
    hold off
    set(gca,'Xlim',[ESDmin ESDmax])
    set(gca,'Ylim',[0.00 1.00])
    xlabel('log (size)')
    ylabel('f (x)')
    grid on
    %...................................................................................
    subplot(2,2,4)
    plot(ESD,flogESD,'k-',ESD,flogESD,'b.')
    hold on
    plot(ESDplus,0,'r*',ESDless,0,'r*')
    hold off
    set(gca,'Xlim',[exp(xmin) exp(xmax)])
    set(gca,'Ylim',[0.00 1.00])
    xlabel('(size)')
    ylabel('f (ESD)')
    grid on
    % sst0 = 15.0; % Must be specified even if mode is ESD (Le Gland, 13/06/2019, to be improved)
%...................................................................................
%===================================================================================
%elseif strcmp(keyTraitAxis,'SST')
%===================================================================================
%...................................................................................
    % yaxis = [ymin:ydel:ymax]; %SST optimal values for each phy species-j
    ymean = mean(yaxis); %SST of the environment.
    %................................................................................... 
    %...................................................................................
    %%fybis = (1.0) * exp(-(yaxis - ymean).^2 / (2*ygamma^2)); %SST limitation [%]
    %...................................................................................
    fy = (1.0) * exp(-(yaxis - ymean).^teta / (2*ygamma^teta)); %Pseudo-Gaussian curve reaching max.value = 1.0 (ie. area bigger than 1.0)
    %...................................................................................
    i = 0;
    sstj = yaxis; 
    sst0 = 15.0;
    Qsst = ones(3,length(fy));
    % $$$ QsstEppley = Q10.^((sstj-sst0)/10); %Temperature limitation at all sst values.
    for ssti = [5:10:25]
        i = i + 1;
        ify = (1.0) * exp(-(yaxis - ssti).^teta / (2*ygamma^teta)); %Pseudo-Gaussian curve reaching max.value = 1.0 (ie. area bigger than 1.0)
        Qsst(i,:) = ify * (Q10a^((ssti - sst0)/10)); %USAR ESTA!!!!
    % $$$     Qsst(i,:) = ify * (Q10^((ssti)/10));
    % $$$     Qsst(i,:) = ify .* QsstEppley;
    end
    %...................................................................................
    figure(100)
    % jamstecrest_gaussecomodel1D_samplingplots(phy_cont,phy_disc2,xrng,yrng,ndays,fignum);
    % $$$ subplot(2,2,1)
    % $$$ plot(yaxis,fy)
    % $$$ grid on
    % $$$ subplot(2,2,2)
    plot(yaxis,Qsst)
    set(gca,'Xlim',[0 30])
    xlabel('y')
    ylabel('Qsst')
    legend('SST = 5','SST = 15','SST = 25')
    grid on
    %...................................................................................
    %%print('-dpng','-r300','jamstecrest_gaussecomodel1D_fig000.png')
%===================================================================================
%end %endif keyTraitAxis == ESD

% Use functions to plot trade-offs (Le Gland, 25/10/2019)
% jamstecrest_uptaketradeoff(amup0,alp0,amup,aalp,[0.03, 0.1, 0.3, 1.0, 3.0],'plot')
% jamstecrest_temptradeoff([20, 24, 28],'plot')
% jamstecrest_tradeoff(1.5,0.5,[0.05 0.2 0.8 3.2],20,1.8,[15:5:30],15);
fignum = 15;
%jamstecrest_tradeoff(mup0,amup,[0.05 0.2 0.8 3.2],temp0,Q10a,[15:5:30],fignum);
%jamstecrest_tradeoff(mup0,amup,[exp(-2.5:1:1.5)],temp0,Q10a,[18:3:30],fignum);
%jamstecrest_tradeoff(mup0,amup,exp(-1.5:0.5:0.5),temp0,Q10a,18:4:30,fignum);
%jamstecrest_tradeoff(mup0,amup,[0.1,0.5,2.0],temp0,Q10a,18:4:30,fignum);
%jamstecrest_tradeoff(mup0,amup,[0.1,0.5,2.0],temp0,Q10a,18:4:30,fignum);
%return


%%%%%%%%%%%%%%%%%%%%%%%
%PLOT UPTAKE TRADE OFF:
%%%%%%%%%%%%%%%%%%%%%%%
% %if strcmp(keyTraitAxis,'ESD')
%     %===================================================================================
%     %...................................................................................
%     %%DINaxis = 10.^[-3:1:1]; %[mmolN*m-3] (for ESDmin = 0.02 and ESDmax = 200)
%     %...................................................................................
%     DINaxis = 10.^[-1:1:2]; %[mmolN*m-3]
%     %...................................................................................
%     [junkarg,DINaxis] = jamstecrest_uptaketradeoff_old(log([0.1,0.5,2,10,50]),'Lanimal');
%     %...................................................................................
%     LXaxis  = ones(length(DINaxis),length(xaxis))*nan;
%     QXaxis  = ones(length(DINaxis),length(xaxis))*nan;
%     UXaxis  = ones(length(DINaxis),length(xaxis))*nan;
%     UXaxisMax  = ones(length(DINaxis),1)*nan;
%     %...................................................................................
%     % $$$ alp0 = (mup0 / knp0); 
%     % $$$ aal = (amup - aknp);
%     %...................................................................................
%     % $$$ MUX = mup0 * exp(amup*xaxis);
%     % $$$ KSX = knp0 * exp(aknp*xaxis);
%     % $$$ ALX = alp0 * exp(aal*xaxis);
%     %...................................................................................
%     % $$$ for i = 1:length(DINaxis)
%     % $$$ %    UXbis(i,:) = MUX .* ( DINaxis(i) ./ (KSX + DINaxis(i)) );
%     % $$$     UXbis(i,:) = ALX .* ( DINaxis(i) ./ (1 + DINaxis(i)./KSX) );
%     % $$$ end
%     %...................................................................................
%     %%[MUX,KSX] = jamstecrest_uptaketradeoff(xaxis,'Vallina');
%     [MUX,KSX] = jamstecrest_uptaketradeoff_old(xaxis,'Lanimal');
%     %...................................................................................
%     showESDMUX = [ESD(:),MUX(:)]
%     %...................................................................................
%     for i = 1:length(DINaxis)
%         LXaxis(i,:) = (KSX ./ (KSX + DINaxis(i)));
%         QXaxis(i,:) = (DINaxis(i) ./ (DINaxis(i) + KSX));
%         UXaxis(i,:) = MUX .* QXaxis(i,:); 
%         UXaxisMax(i,1) = max(UXaxis(i,:));
%     end
%     UXaxisStar = UXaxis ./ (UXaxisMax*ones(1,length(xaxis)));
%     %...................................................................................
%     myXtick = [0.02,0.10,0.50,2.00,10,50,200];
%     %...................................................................................
%     figure(15)
%     %...................................................................................
%     subplot(2,2,1)
%     plot(xaxis,UXaxis(1,:),'g-')
%     hold on
%     plot(xaxis,UXaxis(2,:),'r-')
%     hold on
%     plot(xaxis,UXaxis(3,:),'k-')
%     hold on
%     plot(xaxis,UXaxis(4,:),'b-')
%     hold on
%     plot(xaxis,UXaxis(5,:),'c-')
%     hold off
%     legend([['DIN = ';'DIN = ';'DIN = ';'DIN = ';'DIN = '],num2str(DINaxis(:))],'Location','NorthWest');
%     set(gca,'Ylim',[0 4.0])
%     set(gca,'Xlim',[xmin xmax],'Xtick',log(myXtick),'Xticklabel',(myXtick))
%     xlabel('ESD [um]')
%     ylabel('Growth rate [d-1]')
%     title('Uptake')
%     grid on
%     %...................................................................................
%     subplot(2,2,2)
%     plot(xaxis,UXaxisStar(1,:),'g-')
%     hold on
%     plot(xaxis,UXaxisStar(2,:),'r-')
%     hold on
%     plot(xaxis,UXaxisStar(3,:),'k-')
%     hold on
%     plot(xaxis,UXaxisStar(4,:),'b-')
%     hold on
%     plot(xaxis,UXaxisStar(5,:),'c-')
%     hold off
%     set(gca,'Ylim',[0 1.0])
%     set(gca,'Xlim',[xmin xmax],'Xtick',log(myXtick),'Xticklabel',(myXtick))
%     xlabel('ESD [um]')
%     ylabel('Growth factor [n.d.]')
%     title('Uptake normalized')
%     grid on
%     %...................................................................................
%     subplot(2,2,3)
%     semilogy(xaxis,UXaxis(1,:),'g-')
%     hold on
%     semilogy(xaxis,UXaxis(2,:),'r-')
%     hold on
%     semilogy(xaxis,UXaxis(3,:),'k-')
%     hold on
%     semilogy(xaxis,UXaxis(4,:),'b-')
%     hold on
%     semilogy(xaxis,UXaxis(5,:),'c-')
%     hold off
%     %%set(gca,'Ylim',[0 5.0])
%     %%set(gca,'Ylim',[0 3.0])
%     set(gca,'Xlim',[xmin xmax],'Xtick',log(myXtick),'Xticklabel',(myXtick))
%     xlabel('ESD [um]')
%     ylabel('Growth rate [d-1]')
%     title('Uptake')
%     grid on
%     %...................................................................................
%     subplot(2,2,4)
%     plot(exp(xaxis),UXaxisStar(1,:),'g-')
%     hold on
%     plot(exp(xaxis),UXaxisStar(2,:),'r-')
%     hold on
%     plot(exp(xaxis),UXaxisStar(3,:),'k-')
%     hold on
%     plot(exp(xaxis),UXaxisStar(4,:),'b-')
%     hold on
%     plot(exp(xaxis),UXaxisStar(5,:),'c-')
%     hold off
%     set(gca,'Ylim',[0 1.0])
%     set(gca,'Xlim',[exp(xmin) exp(xmax)])%,'Xtick',(myXtick),'Xticklabel',(myXtick))
%     xlabel('exp(ESD) [um]')
%     ylabel('Growth factor [n.d.]')
%     title('Uptake normalized')
%     grid on
%     %...................................................................................
%     %%print('-dpng','-r300','jamstecrest_gaussecomodel1D_fig000.png')
%     %%return
%     %===================================================================================
%     %...................................................................................
%     ESDvec = [0.5, 1.0, 2.0, 5.0, 10, 20, 50, 100, 200]; 
%     %...................................................................................
%     %%xvec = log( 2*(10.^[-2:1:2]) ); %log([um])
%     %%xvec = log( [0.2, 0.5, 1.0, 2.0, 5.0, 10, 20, 50] ); %log([um])
%     xvec = log(ESDvec); %log([um])
%     %...................................................................................
%     DIN = [0:0.001:2.0]; 
%     %...................................................................................
%     QDIN = ones(length(xvec),length(DIN))*nan;
%     UDIN = ones(length(xvec),length(DIN))*nan;
%     %...................................................................................
%     [MUXvec,KSXvec] = jamstecrest_uptaketradeoff_old(xvec,'Lanimal');
%     %...................................................................................
%     for j = 1:length(xvec)
%         QDIN(j,:) = DIN ./ (KSXvec(j) + DIN);
%         UDIN(j,:) = MUXvec(j) .* QDIN(j,:);
%     end
%     %...................................................................................
%     figure(16)
%     % jamstecrest_gaussecomodel1D_samplingplots(phy_cont,phy_disc2,xrng,yrng,ndays,fignum);
%     plot(DIN,UDIN)
%     grid on
%     pause(1)
%     return
%     close all
%     %...................................................................................
%     %===================================================================================
% %end %endif keyTraitAxis == ESD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIAL CONDITIONS FOR THE CONTINOUS TRAIT MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%if strcmp(keyTraitAxis,'ESD')
    %...............................................................................
    xave0  = log(ESDave); %Phy diameter average [log(um)]
    xstd0  = xsigmaPcnt*(log(ESDmax) - log(ESDmin));
    % xvar0 = (xstd0)^2; %Phy diameter variance [log(um)^2]
    xxvar0 = (xstd0)^2; % Le Gland, 04/07/2019
    %...............................................................................
%elseif strcmp(keyTraitAxis,'SST')
    %...............................................................................
% $$$     xave0 = 5.0; %Phy Topt average [Celsius]
% $$$     xstd0 = 2.0; %Phy Topt deviate [Celsius]
% $$$     xvar0 = xstd0^2; %Phy Topt varianz [Celsius^2]
    %...............................................................................
    %%xave0 = mean(sst); %Phy Topt average [Celsius]
    %xave0 = mean(sst(1:ndays));
    %xstd0 =  2.0; %Phy Topt deviate [Celsius]
    %%xstd0 = 0.002; % Start near zero variance (Le Gland, 12/06/2019)
    %xvar0 = (xstd0)^2; %Phy Topt varianz [Celsius^2]
    % Temperature is now a y trait (Le Gland, 01/07/2019)
    yave0  = mean(sst(1:ndays));
    %ystd0  = 2.0; % Phy Topt deviate [Celsius]
    ystd0 = ysigmaPcnt*(ymax - ymin);
    yyvar0 = (ystd0)^2;
    % Correlation between traits is initialized at zero 
    xycov0 = 0;
    %...............................................................................
%end
%===================================================================================
if strcmp(keyPhysics,'yes')
    %...............................................................................
    xave_star0  = xave0*phy0;
    %xvar_star0 = xvar0*phy0;
    xxvar_star0 = (xxvar0+xave0.^2)*phy0; % Le Gland (05/06/2019)
    %xstd_star0 = xstd0*phy0;
    %...............................................................................
    yave_star0  = yave0*phy0; % Initial raw moments of temperature trait and correlation (Le Gland, 01/07/2019)
    yyvar_star0 = (yyvar0+yave0.^2)*phy0; 
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
BOX0 = box0*vdepths; %column vector with I.C. for Det (ndepths vertical nodes).
%........................................................................
gbeta = gbeta0;
galfa = galfa0*vdepths; %Constant grazing-alfa KTW parameter.
%numut = numut0*vdepths; %Trait diffusion area [trait]^2
% One trait diffusion value for each trait (Le Gland, 10/09/2019)
numutx = numutX0*vdepths;
numuty = numutY0*vdepths;
%........................................................................
if strcmp(keyPhysics,'not')
    %%keyScaling = 'linear';
    %...................................................................................
    keyScaling = 'logarithmic';
    %...................................................................................
    % yexp is ambiguous, I change it to zexp. It is not used later anyway ! (Le Gland, 01/07/2019)
    %[yexp,fexp,fexpesc] = jamstecrest_gaussecomodel1D_fexpfunction(ndepths,keyScaling);
    [zexp,fexp,fexpesc] = jamstecrest_gaussecomodel1D_fexpfunction(ndepths,keyScaling);
    %...................................................................................
    %===================================================================================
    %...................................................................................
    galfa = fexpesc + 1.0;
    %numut = fexpesc * numut0; %Trait diffusion area [trait]^2
    numutx = fexpesc * numutX0;
    numuty = fexpesc * numutY0;
    %........................................................................
    if strcmp(keyGalfaConstant,'yes')
        galfa = galfa0*vdepths; %Constant grazing-alfa KTW parameter.
    end
    %........................................................................
    if strcmp(keyTraitDiffusion,'yes')
        % numut = numut0*vdepths;
        numutx = numutX0*vdepths;
        numuty = numutY0*vdepths;
    end
    %........................................................................
    if strcmp(keyTraitDiffusionConstant,'yes')
        % numut = numut0*vdepths;
        numutx = numutX0*vdepths;
        numuty = numutY0*vdepths;
    end
    %........................................................................
end %endif keyPhysics 
%========================================================================
%........................................................................
%COLUMN VECTOR OF INITIAL CONDITIONS:
%........................................................................
if strcmp(keyPhysics,'not')
    %....................................................................
    %PDF = [PDFave0,PDFvar0,PDFstd0]; %Statistics of Phy size-distribution.
    PDF = [PDFxave0,PDFyave0,PDFxxvar0,PDFyyvar0,PDFxycov0]; % 2 traits (Le Gland, 01/07/2019)
    PDF_K = [PDFxave0,PDFxxvar0]; % Le Gland, 26/11/2019
    PDF_T = [PDFyave0,PDFyyvar0]; % Le Gland, 26/11/2019
    %....................................................................
elseif strcmp(keyPhysics,'yes')
    %....................................................................
    %PDF = [PDFave_star0,PDFvar_star0,PDFstd_star0]; %Statistics of Phy size-distribution.
    PDF = [PDFxave_star0,PDFyave_star0,PDFxxvar_star0,PDFyyvar_star0,PDFxycov_star0]; % 2 traits (Le Gland, 01/07/2019)
    PDF_K = [PDFxave_star0,PDFxxvar_star0]; % Le Gland, 26/11/2019
    PDF_T = [PDFyave_star0,PDFyyvar_star0]; % Le Gland, 26/11/2019
    %....................................................................
end
%........................................................................
BIO = [PHY0,ZOO0,DIN0,PON0,BOX0]; %Biomasses of plankton and nutrients.
%........................................................................
V0 = []; %Initalize vector.
%........................................................................
PDF_V0 = PDF(:); % Transform a matrix (ndepths, 5) into a vector (ndepths * 5)
PDF_V0_K = PDF_K(:); % Le Gland, 26/11/2019
PDF_V0_T = PDF_T(:); % Le Gland, 26/11/2019
%........................................................................
BIO_V0 = BIO(:);
%........................................................................
%V0 = [PDF_V0;BIO_V0]; %column vector with all intitial conditions.
V0 = [BIO_V0;PDF_V0]; % Le Gland, 25/11/2019
V0_K = [BIO_V0;PDF_V0_K]; % Le Gland, 26/11/2019
V0_T = [BIO_V0;PDF_V0_T]; % Le Gland, 26/11/2019
%........................................................................
ntot0 = sum(BIO_V0(:)); %initial total mass over the column water [mmolN*m-3]
%........................................................................
ntot0pernode = ntot0/ndepths; %initial total mass at each grid cell [mmolN*m-3]
%........................................................................
%========================================================================
%%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIAL CONDITIONS FOR THE DISCRETE TRAIT MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------------
%NOTE: Size Spectra Model (SSM) --> Discrete trait model
%------------------------------------------------------------------------
%========================================================================
%........................................................................
%%nphy = 64; %Number of phy species (as many as you wish) (Standard)
%nphy = 61;
%%nphy = 32;
%%nphy = 256; %Number of phy species (as many as you wish)
% Case with 2 traits (Le Gland, 04/07/2019)
%nxphy = 13;
%nyphy = 11;
%........................................................................
nzoo =  1; %Number of zoo species (must be always only one!)
ndin =  1; %DIN (must be always only one!)
npon =  1; %PON (must be always only one!)
nbox =  1; %Virtual box (must be always only one!)
%........................................................................
%========================================================================
%DISCRETE TRAITS:
%if strcmp(keyTraitAxis,'ESD') % Le Gland (26/04/2019)
    %........................................................................
    minESDphy = ESDmin; %0.02 um 
    maxESDphy = ESDmax; %200 um 
    %........................................................................
    % Is base2 mode really needed ? (Le Gland, 29/08/2019)
    % Is sizemin really required if xmin is already present ? (Le Gland, 18/10/2019)
    xsizemin = log(minESDphy); 
    xsizemax = log(maxESDphy); 
    %........................................................................
    %%xsizedel = (xsizemax - xsizemin) / (nphy); 
    %%xsizerng = [xsizemin+xsizedel:xsizedel:xsizemax];
    %........................................................................
    %========================================================================
    %ESD AXIS USING THE VALUES AT THE BORDERS: 
    %........................................................................
    % xsizedel = (xsizemax - xsizemin) / (nphy - 1); %my approach. (used again instead of the former by Le Gland since 04/06/2019)
    xsizedel = (xsizemax - xsizemin) / (nxphy - 1);
    xsizerng = [xsizemin:xsizedel:xsizemax];
    % xtraitdel = (xsizemax - xsizemin) / (nphy - 1); %my approach. (used again instead of the former by Le Gland since 04/06/2019)
    % xtraitrng = [xsizemin:xtraitdel:xsizemax];
    % $$$ %........................................................................
    % $$$ % $$$ ESDlogspace10 = logspace(log10(ESDmin),log10(ESDmax),nphy); %Lan Smith approach.
    % $$$ % $$$ xsizerng = log10(ESDlogspace10) ./ log10(exp(1)); 
    %........................................................................
    %========================================================================
    %ESD AXIS USING THE VALUES AT THE MIDDLE POINT BETWEEN BORDERS: 
    %........................................................................
    %%xsizedel = (xsizemax - xsizemin) / (nphy); %my approach.
    %%xsizerng = [xsizemin+(0.5*xsizedel):xsizedel:xsizemax-(0.5*xsizedel)];
    %........................................................................
    % $$$ ESDlogspace10 = logspace(log10(ESDmin),log10(ESDmax),(nphy+1)); %Lan Smith approach.
    % $$$ xsizeseg = log10(ESDlogspace10) ./ log10(exp(1)); %Segment borders (nphy+1 points) 
    % $$$ xsizerng = xsizeseg(1:end-1) + (1/2)*diff(xsizeseg); %Middel points (nphy points)
    % $$$ xsizesegdiff = diff(xsizeseg); %All segments should be equally spaced. 
    % $$$ xsizedel = xsizesegdiff(1); %So I can use the spacing of just first the segment. 
    %........................................................................
    %========================================================================
    %........................................................................
    logESDphy = xsizerng;
    %logESDphy = xtraitrng; % Le Gland, 04/06/2019
    %........................................................................
    ESDphy = exp(logESDphy); 
    %........................................................................
    %xtrait = logESDphy(1:nphy);
    xtrait = logESDphy(1:nxphy); % Le Gland (04/07/2019)
%........................................................................
%========================================================================
%........................................................................
%elseif strcmp(keyTraitAxis,'SST')
    %%ytoptmin = -2.0; %[Celsius]
    %%ytoptmax = 32.0; %[Celsius]
    % Keeping optimal temperatures between sstmin and sstmax makes more sense (Le Gland, 04/06/2019)
    % ytoptmin = sstmin;
    % ytoptmin = 0;
    ytoptmin = ymin;
    % ytoptmax = sstmax;
    % ytoptmax = 30;
    ytoptmax = ymax;
    ytoptdel = (ytoptmax - ytoptmin) / (nyphy - 1); 
    %........................................................................
    TOPT = [ytoptmin:ytoptdel:ytoptmax];
    %........................................................................
    ytrait = TOPT(1:nyphy);
    %ytrait = TOPT(1:nphy); % Called ytrait again since 01/07/2019 (Le Gland)
    %xtrait = TOPT(1:nphy); % Calling it xtrait is required for following steps (Le Gland, 04/06/2019)
    %........................................................................
    %xmean = (xmin + xmax) / 2;
    %xsigma = 0.2*(xmax-xmin); % 0.2 for consistency with the continuous model
    %%xsigma = 0.0002*(xmax-xmin); % Start at zero variance (Le Gland, 12/06/2019)
    %xplus = xmean + xsigma;
    %xless = xmean - xsigma;
    ymean = (ymin + ymax) / 2;
    % ysigma = 0.2*(ymax-ymin);
    ysigma = ysigmaPcnt*(ymax-ymin); % Le Gland, 18/09/2019 
    yplus = ymean + ysigma;
    yless = ymean - ysigma;
    %........................................................................
    %xdelta = xplus - xless; 
    ydelta = yplus - yless;
    %........................................................................
    %xtraitdel = (ytoptmax - ytoptmin) / (nphy - 1);
    %xtraitrng = [ytoptmin:xtraitdel:ytoptmax];
    ytraitdel = (ytoptmax - ytoptmin) / (nyphy - 1);
    ytraitrng = [ytoptmin:ytraitdel:ytoptmax];
    %........................................................................
    %========================================================================
    %........................................................................
%end % end if strcmp(keyTraitAxis,'ESD') (Le Gland, 26/04/2019)
% if strcmp(keyTraitAxis,'ESD') % Le Gland (26/04/2019)
% use expressions valid for all traits (Le Gland, 04/06/2019)
% $$$ flogESDphy = (1.0 / (xsigma * sqrt(2*pi))) * exp( -(xtrait - xmean).^2 / (2*xsigma^2) );
%if xsigma == 0;
%    ftraitphy = ( xtrait == xmean ); % If no variance, all plankton is concentrated at xmean (Le Gland, 12/06/2019)
%else
%ftraitphy = (1.0 / (xsigma * sqrt(2*pi))) * exp( -(xtrait - xmean).^2 / (2*xsigma^2) );
%end
% Case with 2 traits (Le Gland, 04/07/2019)
fxtraitphy = (1.0 / (xsigma * sqrt(2*pi))) * exp( -(xtrait - xmean).^2 / (2*xsigma^2) );
fytraitphy = (1.0 / (ysigma * sqrt(2*pi))) * exp( -(ytrait - ymean).^2 / (2*ysigma^2) );
% fxytraitphy is a multivariate normal distribution without correlation (Le Gland, 04/07/2019)    
fxytraitphy = fxtraitphy'*fytraitphy;
%........................................................................
% $$$ xtrait1 = xtrait(1:end-1);
% $$$ xtrait2 = xtrait(2:end);
% $$$ xtraitCentered = (1/2)*(xtrait1 + xtrait2);
% $$$ flogESDphy = (1.0 / (xsigma * sqrt(2*pi))) * exp( -(xtraitCentered - xmean).^2 / (2*xsigma^2) );
%........................................................................
%========================================================================
figure(2020)
%...................................................................................
subplot(2,2,1)
% $$$ plot(xtrait,flogESDphy,'k-',xtrait,flogESDphy,'b.')
% $$$ plot(xtrait,ftraitphy,'k-',xtrait,ftraitphy,'b.')
plot(xtrait,fxtraitphy,'k-',xtrait,fxtraitphy,'b.')
hold on
plot(xless,0,'r*',xplus,0,'r*')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('log (size)')
ylabel('f (x)')
grid on
%...................................................................................
%if strcmp(keyTraitAxis,'ESD')
    subplot(2,2,3)
    % $$$ plot(ESDphy,flogESDphy,'k-',ESDphy,flogESDphy,'b.')
    % $$$ plot(ESDphy,ftraitphy,'k-',ESDphy,ftraitphy,'b.')
    plot(ESDphy,fxtraitphy,'k-',ESDphy,fxtraitphy,'b.')
    hold on
    plot(ESDless,0,'r*',ESDplus,0,'r*')
    hold off
    set(gca,'Xlim',[ESDmin ESDmax])
    set(gca,'Ylim',[0.00 1.00])
    xlabel('log (size)')
    ylabel('f (x)')
    grid on
    %...................................................................................
%end % if strcmp(keyTraitAxis,'ESD') % Le Gland (26/04/2019)
pause(1)
close all

%========================================================================
%SINGLE VALUE:
% This version (same initial concentration for all species) is not used
% (Le Gland, 04/11/2019)
%........................................................................
%phydisc0 = (phy0/nphy); %Phyplankton [mmolN*m-3]
phydisc0 = (phy0/(nxphy*nyphy)); % Case with 2 traits (Le Gland, 04/07/2019)
zoodisc0 = (zoo0/nzoo); %Zooplankton [mmolN*m-3]
dindisc0 = (din0/ndin); %Nutrient [mmolN*m-3]
pondisc0 = (pon0/npon); %Detritus [mmolN*m-3]
boxdisc0 = (box0/nbox); %Virtual box [mmolN*m-3]
%........................................................................
%PHYdisc0 = phydisc0 * ones(ndepths*nphy,1);
PHYdisc0 = phydisc0 * ones(ndepths*nxphy*nyphy,1); % Le Gland, 04/07/2019
ZOOdisc0 = zoodisc0 * ones(ndepths*nzoo,1);
DINdisc0 = dindisc0 * ones(ndepths*ndin,1);
PONdisc0 = pondisc0 * ones(ndepths*npon,1);
BOXdisc0 = boxdisc0 * ones(ndepths*nbox,1);
%........................................................................
%========================================================================
%FOR GAUSSIAN DISTRIBUTION OF INITIAL CONDITIONS OF PHYTOPLANKTON: 
%........................................................................
%%phydisc0gaussian = phy0 * (flogESDphy .* xsizedel); %Phy size distribution at t0. 
%%phydisc0gaussian = phy0 * (ftraitphy .* xtraitdel); % Generic trait (Le Gland, 04/06/2019)
%phydisc0gaussian = phy0 * (fxtraitphy .* xsizedel); % Le Gland, 01/07/2019
%phydisc0gaussian = phy0 * (fxytraitphy .* (xsizedel*ytraitdel));% I try to replace it by a matrix (Le Gland, 04/07/2019)
phydisc0gaussian = max(eps, phy0 * (fxytraitphy .* (xsizedel*ytraitdel))); % Avoid strictly zero values (Le Gland, 11/11/2019)
% This formula is correct only for zero initial correlation
%%phydisc0gaussian = phy0 * (); % 2 traits (Le Gland, 01/07/2019, to see later)

%%phydisc0gaussian = phy0 * (flogESDphy); %I THINK THIS IS THE CORRECT ONE!!! (to be independent of x-trait resolution)
%........................................................................
% PHYdisc0gaussian1D = []; 
% for i = 1:ndepths 
%     PHYdisc0gaussian1D = [PHYdisc0gaussian1D;phydisc0gaussian]; 
% end
% In the case with multiple traits, PHYdisc0gausssian1D has to be a matrix
% (Le Gland, 04/07/2019)
PHYdisc0gaussian1D = zeros(ndepths,nxphy,nyphy);
for i = 1:ndepths 
    PHYdisc0gaussian1D(i,:,:) = phydisc0gaussian;
end

%........................................................................
PHYdisc0 = PHYdisc0gaussian1D(:); % I think it is less confusing to work
% with the matrix, especially in the multivariate case (Le Gland, 04/07/2019)
% But this might be a problem in V0 !
%........................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%test: 
%........................................................................
%[xmeanBis,xsigmaBis] = jamstecrest_geometricmean(phydisc0gaussian,logESDphy); 
%........................................................................
%xmeanTest = [xmean,xmeanBis] %They should be the same. 
%xsigmaTest = [xsigma,xsigmaBis] %They should be the same. 
%........................................................................
figure(10)
%plot(xaxis,fx,'-b.')
plot(xaxis,fxtraitphy,'-b')
hold on
%plot(logESDphy,flogESDphy,'r*')
plot(xtrait,fxtraitphy,'r*')
hold off
grid on
%........................................................................
pause(0.5) 
close all 
pause(1)
%%return
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
%mlength = length(PDF_V0);
mlength = length(BIO_V0); % Le Gland, 25/11/2019
%........................................................................
% $$$ Iave = [1:ndepths];
% $$$ Ivar = [(ndepths*1)+1:(ndepths*2)];
% $$$ Istd = [(ndepths*2)+1:(ndepths*3)];
% For 2 traits (Le Gland, 01/07/2019)
%Ixave  = [1:ndepths];
%Iyave  = [(ndepths*1)+1:(ndepths*2)];
%Ixxvar = [(ndepths*2)+1:(ndepths*3)];
%Iyyvar = [(ndepths*3)+1:(ndepths*4)];
%Ixycov = [(ndepths*4)+1:(ndepths*5)];
%........................................................................
%Iphy = mlength + [1:ndepths];
%Izoo = mlength + [(ndepths*1)+1:(ndepths*2)];
%Idin = mlength + [(ndepths*2)+1:(ndepths*3)];
%Ipon = mlength + [(ndepths*3)+1:(ndepths*4)];
%Ibox = mlength + [(ndepths*4)+1:(ndepths*5)];
%........................................................................
% Change in indices to accomodate both 1 and 2-trait models (Le Gland, 25/11/2019)
Iphy = 1:ndepths;
Izoo = (ndepths*1)+1:(ndepths*2);
Idin = (ndepths*2)+1:(ndepths*3);
Ipon = (ndepths*3)+1:(ndepths*4);
Ibox = (ndepths*4)+1:(ndepths*5);
%........................................................................
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
%wIndex = [Iave;Ivar;Istd;Iphy;Izoo;Idin;Ipon;Ibox]
% For 2 traits (Le Gland, 01/07/2019)
%wIndex = [Ixave;Iyave;Ixxvar;Iyyvar;Ixycov;Iphy;Izoo;Idin;Ipon;Ibox]
% Moments are pushed to the end (Le Gland, 25/11/2019)
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
    Jvarj = [(ndepths*(nsum1))+1:(ndepths*(nsum2))];
    %....................................................................
    Jstr = ['J',varname{j}];
    %....................................................................
% $$$     assignin('base',Jstr,Jvarj); %NO USAR ESTAR FORMA.
    %....................................................................
    myassign = [Jstr,' = Jvarj;']; %USAR ESTA FORMA MEJOR FOR GLOBAL VARIABLES.
    %....................................................................
    eval(myassign)
    %....................................................................
end

%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRAIT DIFFUSION FOR DISCRETE MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
if strcmp(keyTraitDiffusion,'yes')
%...................................................................................
    %dx = xsizedel;
    %dx = xtraitdel; % Le Gland, 04/06/2019
    % Case with 2 traits (Le Gland, 02/07/2019)
    dx = xdel;
    dy = ydel;
    %...................................................................................
    % $$$ muXpercapita = muXmax/ntot0; %Mutation specific rate: [d-1] / [mmolN * m-3]
    % $$$ kppxStarMax = muXmax*sigmaXm^2 ./ ntot0; %Normalized mutation rate: [d-1] x [trait]^2 / [mmolN * m-3] = [trait^2 * day-1] / [mmolN * m-3]
    %...................................................................................
    % $$$ kppxStar = ones(1,nphy) * kppxStarMax; %[trait^2 * day-1] / [mmolN * m-3] %Constant with size
    %...................................................................................
    % $$$ dx = xsizedel; 
    % $$$ xrng  = xtrait; 
    % $$$ xrngI = [xrng(1)+0.5*dx:dx:xrng(end)-0.5*dx];
    % $$$ kppxStarI = interp1(xrng,kppxStar,xrngI);
    %...................................................................................
    %VERTICALLY CONSTANT FOR FULL 1D SIMULATIONS:
    % $$$ KPPXstar  = vdepths*kppxStar ; %[depth,cellsize]
    % $$$ KPPXstarI = vdepths*kppxStarI; %[depth,cellsize]
    %...................................................................................
    %VERTICAL GRADIENT FOR PARALLEL 0D SIMULATIONS:
    % $$$ %%KPPXstar  = (vdepths*kppxStar ) .* (fexpesc*ones(1,nphy)); %[depth,cellsize]
    % $$$ %%KPPXstarI = (vdepths*kppxStarI) .* (fexpesc*ones(1,nphy-1)); %[depth,cellsize]
    %...................................................................................
    % $$$ maxKPP = kppxStarMax*ntot0; %[m2*d-1]
    %...................................................................................
    %===================================================================================
    %...................................................................................
    %if strcmp(keyTraitAxis,'ESD')
    %    mup_mean = mup0 .* exp(amup*mean(xtrait)); %[d-1]
    %elseif strcmp(keyTraitAxis,'SST')
    %    mup_mean = mup0 .* Q10^(mean(xtrait)/10);
    %end
    mup_mean = mup0 .* exp(amup*mean(xtrait)) .* Q10a^(mean(ytrait)/10); % Case with 2 traits (Le Gland, 02/07/2019)
    %rx = mup_mean; %[d-1]
    %nux = numut0; %[trait]^2
    % One diffudivity for each trait (Le Gland, 10/09/2019)
    nux = numutX0;
    nuy = numutY0;
    %kppx = (nux.*rx); %[trait]^2 x [d-1] = [trait^2 * d-1] 
    kppx = nux.*mup_mean;
    kppy = nuy.*mup_mean;
    %maxKPP = max(kppx(:)); %[trait^2 * d-1]
    maxKPPx = max(kppx(:));
    maxKPPy = max(kppy(:));
    %...................................................................................
    %===================================================================================
    %...................................................................................
    %StabilityKPP = min(dx^2./(2*maxKPP))
    %StabilityKPPx = min(dx^2./(2*maxKPPx));
    %StabilityKPPy = min(dy^2./(2*maxKPPy));
    StabilityKPP  = min(dx^2./(2*maxKPPx),dy^2./(2*maxKPPy))
    %...................................................................................
    % This could also be improved by an implicit or split explicit time
    % stepping (Le Gland, 18/10/2019). But is it worth it ? It is not
    % realistic with several traits !
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
%...................................................................................
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
%jjday = 1;
%...................................................................................
tcounter = 0;
jcounter = 0;
%...................................................................................
Voutcont = zeros(length(tspan),length(V0));
%...................................................................................
ode45options = odeset('AbsTol',1e-12,'RelTol',1e-6);
%...................................................................................
toc
tic
[Voutcont] = ode4(@jamstecrest_gaussecomodel1D_ode45eqs,tspan,V0);
toc
%...................................................................................
Voutcont = Voutcont'; %Continous Based Model (continous).
% Outputs from 1-trait models (Le Gland, 26/11/2019)
tic
[Voutcont_K] = ode4(@jamstecrest_gaussecomodel1D_ode45eqs_K,tspan,V0_K);
%Voutcont_K = ones(nsteps,140);
toc
Voutcont_K = Voutcont_K';
tic
[Voutcont_T] = ode4(@jamstecrest_gaussecomodel1D_ode45eqs_T,tspan,V0_T);
%Voutcont_T = ones(nsteps,140);
toc
Voutcont_T = Voutcont_T';
%...................................................................................
%save('jamstecrest_gaussecomodel1D_Voutcont.mat','Voutcont')
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
%jjday = 1;
%...................................................................................
icounter = 1;
%...................................................................................
Voutdisc = zeros(length(tspan),length(Vdisc0));
%...................................................................................
ode45options = odeset('AbsTol',1e-12,'RelTol',1e-6);
%...................................................................................
[Voutdisc] = ode4(@jamstecrest_discretemodel1D_ode45eqs,tspan,Vdisc0);
%Voutdisc = ones(nsteps,80+20*nxphy*nyphy);
%...................................................................................
Voutdisc = Voutdisc'; %Size spectra model (discrete).
%...................................................................................
save('jamstecrest_gaussecomodel1D_Voutdisc.mat','Voutdisc')
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
[Vodecont,todecont] = jamstecrest_dailyAve(Voutcont,deltat);
% 1-trait models (Le Gland, 28/11/2019)
[Vodecont_K,todecont] = jamstecrest_dailyAve(Voutcont_K,deltat);
[Vodecont_T,todecont] = jamstecrest_dailyAve(Voutcont_T,deltat);
%...................................................................................
[SDINode,junkarg] = jamstecrest_dailyAve(Sdin,deltat);
%...................................................................................
%===================================================================================
%...................................................................................
% Useless ? (these are only the derivatives) (Le Gland, 26/11/2019) 
% Xavedot  = jamstecrest_dailyAve(Xavedotout,deltat);
% Case with 2 traits (Le Gland, 17/07/2019)
% Yavedot  = jamstecrest_dailyAve(Yavedotout,deltat);
% XXvardot = jamstecrest_dailyAve(XXvardotout,deltat);
% YYvardot = jamstecrest_dailyAve(YYvardotout,deltat);
% XYcovdot = jamstecrest_dailyAve(XYcovdotout,deltat);
% Values for 1-trait models (Le Gland, 26/11/2019)
%...................................................................................
% Case with 2 traits (Le Gland, 17/07/2019)
UXYode = jamstecrest_dailyAve(UXYout,deltat);
GXYode = jamstecrest_dailyAve(GXYout,deltat);
% Values for 1-trait models (Le Gland, 26/11/2019)
UXode = jamstecrest_dailyAve(UXout,deltat);
UYode = jamstecrest_dailyAve(UYout,deltat);
GXode = jamstecrest_dailyAve(GXout,deltat);
GYode = jamstecrest_dailyAve(GYout,deltat);
%...................................................................................
%d1UXdxode = jamstecrest_dailyAve(d1UXdxout,deltat);
%d1GXdxode = jamstecrest_dailyAve(d1GXdxout,deltat);
% Case with 2 traits (Le Gland, 17/07/2019)
% d1UXYdxode = jamstecrest_dailyAve(d1UXYdxout,deltat);
% d1GXYdxode = jamstecrest_dailyAve(d1GXYdxout,deltat);
% d1UXYdyode = jamstecrest_dailyAve(d1UXYdyout,deltat);
% d1GXYdyode = jamstecrest_dailyAve(d1GXYdyout,deltat);
%...................................................................................
%d2UXdxode = jamstecrest_dailyAve(d2UXdxout,deltat);
%d2GXdxode = jamstecrest_dailyAve(d2GXdxout,deltat);
% Case with 2 traits (Le Gland, 17/07/2019)
% d2UXYdxdxode = jamstecrest_dailyAve(d2UXYdxdxout,deltat);
% d2GXYdxdxode = jamstecrest_dailyAve(d2GXYdxdxout,deltat);
% d2UXYdydyode = jamstecrest_dailyAve(d2UXYdydyout,deltat);
% d2GXYdydyode = jamstecrest_dailyAve(d2GXYdydyout,deltat);
% d2UXYdxdyode = jamstecrest_dailyAve(d2UXYdxdyout,deltat);
% d2GXYdxdyode = jamstecrest_dailyAve(d2GXYdxdyout,deltat);
%...................................................................................
% todedot = jamstecrest_dailyAve(todedotout,deltat);
%...................................................................................
%===================================================================================
%...................................................................................
uphyodedisc = uphydaydisc; %Phy specific uptake rate of discrete model [d-1]
%...................................................................................
gphyodedisc = gphydaydisc; 
%...................................................................................
FPHYodedisc = FPHYdaydisc; %Zoo uptake rate of discrete model [mmolN*m-3*d-1]
%...................................................................................
GPHYodedisc = GPHYdaydisc; 
%...................................................................................
%===================================================================================
% $$$ %...................................................................................
% $$$ [uphyodedisc,junkarg] = jamstecrest_dailyAve(uphyoutdisc,deltat); %Phy specific uptake rate of discrete model [d-1]
% $$$ %...................................................................................
% $$$ [gphyodedisc,junkarg] = jamstecrest_dailyAve(gphyoutdisc,deltat);
% $$$ %...................................................................................
% $$$ [FPHYodedisc,junkarg] = jamstecrest_dailyAve(FPHYoutdisc,deltat); %Zoo uptake rate of discrete model [mmolN*m-3*d-1]
% $$$ %...................................................................................
% $$$ [GPHYodedisc,junkarg] = jamstecrest_dailyAve(GPHYoutdisc,deltat);
% $$$ %...................................................................................
%===================================================================================
% %...................................................................................
% [FPHYTodecont,junkarg] = jamstecrest_dailyAve(FPHYToutcont,deltat);
% %...................................................................................
% [EPHYTodecont,junkarg] = jamstecrest_dailyAve(EPHYToutcont,deltat);
% %...................................................................................
% [GPHYTodecont,junkarg] = jamstecrest_dailyAve(GPHYToutcont,deltat);
% %...................................................................................
% [MPHYTodecont,junkarg] = jamstecrest_dailyAve(MPHYToutcont,deltat);
% %...................................................................................
% [FZOOodecont,junkarg] = jamstecrest_dailyAve(FZOOoutcont,deltat);
% %...................................................................................
% [EZOOodecont,junkarg] = jamstecrest_dailyAve(EZOOoutcont,deltat);
% %...................................................................................
% [MZOOodecont,junkarg] = jamstecrest_dailyAve(MZOOoutcont,deltat);
% %...................................................................................
% [FDINodecont,junkarg] = jamstecrest_dailyAve(FDINoutcont,deltat);
% %...................................................................................
% [FPONodecont,junkarg] = jamstecrest_dailyAve(FPONoutcont,deltat);
% %...................................................................................
% %===================================================================================
% %...................................................................................
% [FPHYTodedisc,junkarg] = jamstecrest_dailyAve(FPHYToutdisc,deltat);
% %...................................................................................
% [EPHYTodedisc,junkarg] = jamstecrest_dailyAve(EPHYToutdisc,deltat);
% %...................................................................................
% [GPHYTodedisc,junkarg] = jamstecrest_dailyAve(GPHYToutdisc,deltat);
% %...................................................................................
% [MPHYTodedisc,junkarg] = jamstecrest_dailyAve(MPHYToutdisc,deltat);
% %...................................................................................
% [FZOOodedisc,junkarg] = jamstecrest_dailyAve(FZOOoutdisc,deltat);
% %...................................................................................
% [EZOOodedisc,junkarg] = jamstecrest_dailyAve(EZOOoutdisc,deltat);
% %...................................................................................
% [MZOOodedisc,junkarg] = jamstecrest_dailyAve(MZOOoutdisc,deltat);
% %...................................................................................
% [FDINodedisc,junkarg] = jamstecrest_dailyAve(FDINoutdisc,deltat);
% %...................................................................................
% [FPONodedisc,junkarg] = jamstecrest_dailyAve(FPONoutdisc,deltat);
% %...................................................................................
%===================================================================================
%...................................................................................
%%return
deltaday = 1;
%...................................................................................
%%Jstep = [1:nsteps]; %For whole simulation (all years at all time steps)
%...................................................................................
Jstep = [1:(ndays/deltat)]; %For first year only (at all time steps) 
%...................................................................................
%%Jstep = [(ndays/deltat)*(nyear-1)+1:(ndays/deltat)*nyear]; %For last year only (at all time steps) 
%...................................................................................
Jdays = [(ndays/deltaday)*(nyear-1)+1:(ndays/deltaday)*nyear]; %Last year days (ie. from day 721 until day 1080)
%...................................................................................
%===================================================================================
%...................................................................................
% PHYThdpcont = Voutcont(Iphy,Jstep); %HDP = High Density Points.
% ZOOhdpcont = Voutcont(Izoo,Jstep);
% DINhdpcont = Voutcont(Idin,Jstep);
% PONhdpcont = Voutcont(Ipon,Jstep);
% % 1-trait models (Le Gland, 28/11/2019)
% PHYThdpcont_K = Voutcont_K(Iphy,Jstep);
% ZOOhdpcont_K  = Voutcont_K(Izoo,Jstep);
% DINhdpcont_K  = Voutcont_K(Idin,Jstep);
% PONhdpcont_K  = Voutcont_K(Ipon,Jstep);
% PHYThdpcont_T = Voutcont_T(Iphy,Jstep);
% ZOOhdpcont_T  = Voutcont_T(Izoo,Jstep);
% DINhdpcont_T  = Voutcont_T(Idin,Jstep);
% PONhdpcont_T  = Voutcont_T(Ipon,Jstep);
%...................................................................................
% NTOThdpcont = PHYThdpcont + ZOOhdpcont + DINhdpcont + PONhdpcont;
%...................................................................................
%===================================================================================
if strcmp(keyPhysics,'not')
    %...............................................................................
    %XAVE_hdpcont = Voutcont(Iave,Jstep);
    %XVAR_hdpcont = Voutcont(Ivar,Jstep);
    %XSTD_hdpcont_bis = Voutcont(Istd,Jstep);
    % Case with 2 traits (Le Gland, 17/07/2019)
    % XAVE_hdpcont  = Voutcont(Ixave,Jstep);
    % YAVE_hdpcont  = Voutcont(Iyave,Jstep);
    % XXVAR_hdpcont = Voutcont(Ixxvar,Jstep);
    % YYVAR_hdpcont = Voutcont(Iyyvar,Jstep);
    % XYCOV_hdpcont = Voutcont(Ixycov,Jstep);
    % 1-trait models (Le Gland, 28/11/2019)
    % XAVE_hdpcont_K  = Voutcont_K(Ixave_K,Jstep);
    % XXVAR_hdpcont_K = Voutcont_K(Ixxvar_K,Jstep);
    % YAVE_hdpcont_T  = Voutcont_T(Iyave_T,Jstep);
    % YYVAR_hdpcont_T = Voutcont_T(Iyyvar_T,Jstep);
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    %XAVE_STAR_hdpcont = Voutcont(Iave,Jstep);
    %XVAR_STAR_hdpcont = Voutcont(Ivar,Jstep);
    %XSTD_STAR_hdpcont = Voutcont(Istd,Jstep);
    % Case with 2 traits (Le Gland, 17/07/2019)
    % XAVE_STAR_hdpcont  = Voutcont(Ixave,Jstep);
    % YAVE_STAR_hdpcont  = Voutcont(Iyave,Jstep);
    % XXVAR_STAR_hdpcont = Voutcont(Ixxvar,Jstep);
    % YYVAR_STAR_hdpcont = Voutcont(Iyyvar,Jstep);
    % XYCOV_STAR_hdpcont = Voutcont(Ixycov,Jstep);
    % 1-trait models (Le Gland, 28/11/2019)
    % XAVE_STAR_hdpcont_K  = Voutcont_K(Ixave_K,Jstep);
    % XXVAR_STAR_hdpcont_K = Voutcont_K(Ixxvar_K,Jstep);
    % YAVE_STAR_hdpcont_T  = Voutcont_T(Iyave_T,Jstep);
    % YYVAR_STAR_hdpcont_T = Voutcont_T(Iyyvar_T,Jstep);
    %...............................................................................
    %XAVE_hdpcont = (XAVE_STAR_hdpcont./PHYThdpcont); %mean size
    %XVAR_hdpcont = (XVAR_STAR_hdpcont./PHYThdpcont) - XAVE_hdpcont.^2; %variance
    %XSTD_hdpcont_bis = (XSTD_STAR_hdpcont./PHYThdpcont) - XAVE_hdpcont; %deviation
    % Case with 2 traits (Le Gland, 17/07/2019)
    % XAVE_hdpcont  = (XAVE_STAR_hdpcont./PHYThdpcont);
    % YAVE_hdpcont  = (YAVE_STAR_hdpcont./PHYThdpcont);
    % XXVAR_hdpcont = (XXVAR_STAR_hdpcont./PHYThdpcont) - XAVE_hdpcont.^2;
    % YYVAR_hdpcont = (YYVAR_STAR_hdpcont./PHYThdpcont) - YAVE_hdpcont.^2;
    % XYCOV_hdpcont = (XYCOV_STAR_hdpcont./PHYThdpcont) - XAVE_hdpcont.*YAVE_hdpcont;
    % 1-trait models (Le Gland, 28/11/2019)
    % XAVE_hdpcont_K  = (XAVE_STAR_hdpcont_K./PHYThdpcont_K);
    % XXVAR_hdpcont_K = (XXVAR_STAR_hdpcont_K./PHYThdpcont_K) - XAVE_hdpcont_K.^2;
    % YAVE_hdpcont_T  = (YAVE_STAR_hdpcont_T./PHYThdpcont_T);
    % YYVAR_hdpcont_T = (YYVAR_STAR_hdpcont_T./PHYThdpcont_T) - YAVE_hdpcont_T.^2;
    %...............................................................................
end
%===================================================================================
%...................................................................................
PHYTodecont = Vodecont(Iphy,:);
ZOOodecont = Vodecont(Izoo,:);
DINodecont = Vodecont(Idin,:);
PONodecont = Vodecont(Ipon,:);
% 1-trait models (Le Gland, 28/11/2019)
PHYTodecont_K = Vodecont_K(Iphy,:);
ZOOodecont_K  = Vodecont_K(Izoo,:);
DINodecont_K  = Vodecont_K(Idin,:);
PONodecont_K  = Vodecont_K(Ipon,:);
PHYTodecont_T = Vodecont_T(Iphy,:);
ZOOodecont_T  = Vodecont_T(Izoo,:);
DINodecont_T  = Vodecont_T(Idin,:);
PONodecont_T  = Vodecont_T(Ipon,:);
%...................................................................................
NTOTodecont = PHYTodecont + ZOOodecont + DINodecont + PONodecont;
%...................................................................................
%===================================================================================
if strcmp(keyPhysics,'not')
    %...............................................................................
    %XAVE_odecont = Vodecont(Iave,:);
    %XVAR_odecont = Vodecont(Ivar,:);
    %XSTD_odecont_bis = Vodecont(Istd,:);
    % Case with 2 traits (Le Gland, 17/07/2019)
    XAVE_odecont  = Vodecont(Ixave,:);
    YAVE_odecont  = Vodecont(Iyave,:);
    XXVAR_odecont = Vodecont(Ixxvar,:);
    YYVAR_odecont = Vodecont(Iyyvar,:);
    XYCOV_odecont = Vodecont(Ixycov,:);
    % 1-trait models (Le Gland, 28/11/2019)
    XAVE_odecont_K  = Vodecont_K(Ixave_K,:);
    XXVAR_odecont_K = Vodecont_K(Ixxvar_K,:);
    YAVE_odecont_T  = Vodecont_T(Iyave_T,:);
    YYVAR_odecont_T = Vodecont_T(Iyyvar_T,:);
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    %XAVE_STAR_odecont = Vodecont(Iave,:);
    %XVAR_STAR_odecont = Vodecont(Ivar,:);
    %XSTD_STAR_odecont = Vodecont(Istd,:);
    % Case with 2 traits (Le Gland, 17/07/2019)
    XAVE_STAR_odecont  = Vodecont(Ixave,:);
    YAVE_STAR_odecont  = Vodecont(Iyave,:);
    XXVAR_STAR_odecont = Vodecont(Ixxvar,:);
    YYVAR_STAR_odecont = Vodecont(Iyyvar,:);
    XYCOV_STAR_odecont = Vodecont(Ixycov,:);
    % 1-trait models (Le Gland, 28/11/2019)
    XAVE_STAR_odecont_K  = Vodecont_K(Ixave_K,:);
    XXVAR_STAR_odecont_K = Vodecont_K(Ixxvar_K,:);
    YAVE_STAR_odecont_T  = Vodecont_T(Iyave_T,:);
    YYVAR_STAR_odecont_T = Vodecont_T(Iyyvar_T,:);
    %...............................................................................
    %XAVE_odecont = (XAVE_STAR_odecont./PHYTodecont); %mean trait
    %XVAR_odecont = (XVAR_STAR_odecont./PHYTodecont) - XAVE_odecont.^2; %variance
    %XSTD_odecont_bis = (XSTD_STAR_odecont./PHYTodecont) - XAVE_odecont; %deviation
    % Case with 2 traits (Le Gland, 17/07/2019)
    XAVE_odecont  = (XAVE_STAR_odecont./PHYTodecont);
    YAVE_odecont  = (YAVE_STAR_odecont./PHYTodecont);
    XXVAR_odecont = (XXVAR_STAR_odecont./PHYTodecont) - XAVE_odecont.^2;
    YYVAR_odecont = (YYVAR_STAR_odecont./PHYTodecont) - YAVE_odecont.^2;
    XYCOV_odecont = (XYCOV_STAR_odecont./PHYTodecont) - XAVE_odecont.*YAVE_odecont;
    % 1-trait models (Le Gland, 28/11/2019)
    XAVE_odecont_K  = (XAVE_STAR_odecont_K./PHYTodecont_K);
    XXVAR_odecont_K = (XXVAR_STAR_odecont_K./PHYTodecont_K) - XAVE_odecont_K.^2;
    YAVE_odecont_T  = (YAVE_STAR_odecont_T./PHYTodecont_T);
    YYVAR_odecont_T = (YYVAR_STAR_odecont_T./PHYTodecont_T) - YAVE_odecont_T.^2;
    %...............................................................................
end
%===================================================================================
%...................................................................................
PHYTsspcont = PHYTodecont(:,Jdays);
ZOOsspcont = ZOOodecont(:,Jdays);
DINsspcont = DINodecont(:,Jdays);
PONsspcont = PONodecont(:,Jdays);
% 1-trait models (Le Gland, 28/11/2019)
PHYTsspcont_K = PHYTodecont_K(:,Jdays);
ZOOsspcont_K  = ZOOodecont_K(:,Jdays);
DINsspcont_K  = DINodecont_K(:,Jdays);
PONsspcont_K  = PONodecont_K(:,Jdays);
PHYTsspcont_T = PHYTodecont_T(:,Jdays);
ZOOsspcont_T  = ZOOodecont_T(:,Jdays);
DINsspcont_T  = DINodecont_T(:,Jdays);
PONsspcont_T  = PONodecont_T(:,Jdays);
%...................................................................................
NTOTsspcont = PHYTsspcont + ZOOsspcont + DINsspcont + PONsspcont;
%...................................................................................
%===================================================================================
if strcmp(keyPhysics,'not')
    %...............................................................................
    XAVE_sspcont = XAVE_odecont(:,Jdays);
    %XVAR_sspcont = XVAR_odecont(:,Jdays);
    %% XSTD_sspcont_bis = XSTD_odecont(:,Jdays);
    %XSTD_sspcont_bis = XSTD_odecont_bis(:,Jdays);
    % Case with 2 traits (Le Gland, 17/07/2019)
    YAVE_sspcont  = YAVE_odecont(:,Jdays);
    XXVAR_sspcont = XXVAR_odecont(:,Jdays);
    YYVAR_sspcont = YYVAR_odecont(:,Jdays);
    XYCOV_sspcont = XYCOV_odecont(:,Jdays);
    % 1-trait models (Le Gland, 28/11/2019)
    XAVE_sspcont_K  = XAVE_odecont_K(:,Jdays);
    XXVAR_sspcont_K = XXVAR_odecont_K(:,Jdays);
    YAVE_sspcont_T  = YAVE_odecont_T(:,Jdays);
    YYVAR_sspcont_T = YYVAR_odecont_T(:,Jdays);
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    %XAVE_STAR_sspcont = Vodecont(Iave,Jdays);
    %XVAR_STAR_sspcont = Vodecont(Ivar,Jdays);
    %XSTD_STAR_sspcont = Vodecont(Istd,Jdays);
    % Case with 2 traits (Le Gland, 17/07/2019)
    XAVE_STAR_sspcont  = XAVE_STAR_odecont(:,Jdays);
    YAVE_STAR_sspcont  = YAVE_STAR_odecont(:,Jdays);
    XXVAR_STAR_sspcont = XXVAR_STAR_odecont(:,Jdays);
    YYVAR_STAR_sspcont = YYVAR_STAR_odecont(:,Jdays);
    XYCOV_STAR_sspcont = XYCOV_STAR_odecont(:,Jdays);
    % 1-trait models (Le Gland, 28/11/2019)
    XAVE_STAR_sspcont_K  = XAVE_odecont_K(:,Jdays);
    XXVAR_STAR_sspcont_K = XXVAR_odecont_K(:,Jdays);
    YAVE_STAR_sspcont_T  = YAVE_odecont_T(:,Jdays);
    YYVAR_STAR_sspcont_T = YYVAR_odecont_T(:,Jdays);
    %...............................................................................
    %XAVE_sspcont = (XAVE_STAR_sspcont./PHYTsspcont); %mean trait
    %XVAR_sspcont = (XVAR_STAR_sspcont./PHYTsspcont) - XAVE_sspcont.^2; %variance
    %XSTD_sspcont_bis = (XSTD_STAR_sspcont./PHYTsspcont) - XAVE_sspcont; %deviation
    % Case with 2 traits (Le Gland, 17/07/2019)
    % Could be simplified by putting it out of the "if" condition
    XAVE_sspcont = XAVE_odecont(:,Jdays);
    YAVE_sspcont  = YAVE_odecont(:,Jdays);
    XXVAR_sspcont = XXVAR_odecont(:,Jdays);
    YYVAR_sspcont = YYVAR_odecont(:,Jdays);
    XYCOV_sspcont = XYCOV_odecont(:,Jdays);
    % 1-trait models (Le Gland, 28/11/2019)
    XAVE_sspcont_K  = XAVE_odecont_K(:,Jdays);
    XXVAR_sspcont_K = XXVAR_odecont_K(:,Jdays);
    YAVE_sspcont_T  = YAVE_odecont_T(:,Jdays);
    YYVAR_sspcont_T = YYVAR_odecont_T(:,Jdays);
    %...............................................................................
end
%===================================================================================
%OBTAIN STANDARD DEVIATION BY SQUARING THE VARIANCE:
%...................................................................................
% XSTD_hdpcont = (XVAR_hdpcont.^2);
% XSTD_odecont = (XVAR_odecont.^2);
% XSTD_sspcont = (XVAR_sspcont.^2);
% Le Gland (23/04/2019) : STANDARD DEVIATION IS THE SQUARE ROOT OF
% VARIANCE, NOT ITS SQUARE
% XSTD_hdpcont = sqrt(XVAR_hdpcont);
% XSTD_odecont = sqrt(XVAR_odecont);
% XSTD_sspcont = sqrt(XVAR_sspcont);
% Case with 2 traits (Le Gland, 17/07/2019)
% XSTD_hdpcont = sqrt(XXVAR_hdpcont);
XSTD_odecont = sqrt(XXVAR_odecont);
XSTD_sspcont = sqrt(XXVAR_sspcont);
% YSTD_hdpcont = sqrt(YYVAR_hdpcont);
YSTD_odecont = sqrt(YYVAR_odecont);
YSTD_sspcont = sqrt(YYVAR_sspcont);
% 1-trait models (Le Gland, 28/11/2019)
% XSTD_hdpcont_K = sqrt(XXVAR_hdpcont_K);
XSTD_odecont_K = sqrt(XXVAR_odecont_K);
XSTD_sspcont_K = sqrt(XXVAR_sspcont_K);
% YSTD_hdpcont_T = sqrt(YYVAR_hdpcont_T);
YSTD_odecont_T = sqrt(YYVAR_odecont_T);
YSTD_sspcont_T = sqrt(YYVAR_sspcont_T);
% Correlation (covariance normalized by standard deviations, always between -1 and 1) also needs to be assessed
% Correlation is dimensionless (Le Gland, 02/09/2019)
% XYCOR_hdpcont = XYCOV_hdpcont ./ (XSTD_hdpcont .* YSTD_hdpcont);
XYCOR_odecont = XYCOV_odecont ./ (XSTD_odecont .* YSTD_odecont);
XYCOR_sspcont = XYCOV_sspcont ./ (XSTD_sspcont .* YSTD_sspcont);
%...................................................................................
%===================================================================================
%CHANGE SOME VARNAMES: 
%...................................................................................
ESDphysspAveContBis = []; 
ESDphysspStdContBis = []; 
%...................................................................................
% logESDphyhdpAveCont = XAVE_hdpcont; 
% logESDphyhdpStdCont = XSTD_hdpcont; 
%...................................................................................
logESDphysspAveCont = XAVE_sspcont; 
logESDphysspStdCont = XSTD_sspcont; 
% Plot SST (Le Gland, 23/07/2019)
% Why give a new name when the old name is clear enough on what it is ?
% Using the old name is also better when Topt is the trait
% I try to use the old name whenever possible (Le Gland, 30/08/2019) 
% logESDphyhdpAveCont = YAVE_hdpcont; 
% logESDphyhdpStdCont = YSTD_hdpcont; 
%...................................................................................
% logESDphysspAveCont = YAVE_sspcont; 
% logESDphysspStdCont = YSTD_sspcont;
% Give specific name to temperature variables and correlation (Le Gland, 02/09/2019)
% TOPTphyhdpAveCont = YAVE_hdpcont;
% TOPTphyhdpStdCont = YSTD_hdpcont;
%...................................................................................
TOPTphysspAveCont = YAVE_sspcont;
TOPTphysspStdCont = YSTD_sspcont;
%...................................................................................
% Correlations of the continuous model (Le Gland, 03/09/2019)
% phyhdpCorCont = XYCOR_hdpcont;
physspCorCont = XYCOR_sspcont;
% exponential of 2D entropy (Le Gland, 24/10/2019)
% has the dimension of xtrait*ytrait
% "normal" entropy is 1/2 * log(2*pi*e*det(V))
% phyhdpEntCont = 2*pi*exp(1).*XSTD_hdpcont.*YSTD_hdpcont.*sqrt(1-phyhdpCorCont.^2);
physspEntCont = 2*pi*exp(1).*XSTD_sspcont.*YSTD_sspcont.*sqrt(1-physspCorCont.^2);
%...................................................................................
%===================================================================================
%EXACT SOLUTION OF SHANNON ENTROPY FOR THE CONTINOUS LOG-NORMAL DISTRIBUTION: 
%-----------------------------------------------------------------------------------
%NOTE: 
%-----------------------------------------------------------------------------------
% Quintana etal 2008
% H = (1/2) + log(sqrt(2*pi) * sigmax) 
% H = (1/2) + log(sqrt(2*pi) * logxsigma) + logxave 
%-----------------------------------------------------------------------------------
% Entropy lognormal 1D = log(2*pi*logxsigma^2 * exp(logxave + 1/2))
% Entropy lognormal 1D = log(xsigma .* exp(xave + 1/2) * sqrt(2*pi)); % <https://en.wikipedia.org/wiki/Log-normal_distribution> 
%-----------------------------------------------------------------------------------
%...................................................................................
logESDave = logESDphysspAveCont; 
logESDstd = logESDphysspStdCont; 
%...................................................................................
[ESDaveUpper,ESDaveLower,ESDave,ESDstd,ESDcv] = jamstecrest_meansigmaESD(logESDave,logESDstd);
%...................................................................................
ShannonEntropyTheoreticalCont001 = (1/2) + log(sqrt(2*pi) * logESDstd) + logESDave; %WRONG!!!!
ShannonEntropyTheoreticalCont001bis = log(2*pi*logESDstd.^2 .* exp(logESDave + 1/2));
%...................................................................................

%===================================================================================
%EXACT SOLUTION OF SHANNON ENTROPY FOR THE CONTINOUS ABS-NORMAL DISTRIBUTION: 
%-----------------------------------------------------------------------------------
% Entropy absnormal 1D = (1/2) * ln((2*pi)*(sigma^2)*exp(1)) 
% Entropy absnormal 1D = (1/2) * ln((2*pi)*(sigma^2)) + (1/2) 
%-----------------------------------------------------------------------------------
% Entropy absnormal 1D = (1/2) * (1 + ln(2*pi)) + (1/2) * ln(sigmax^2);
% Entropy absnormal 1D = (1/2) * (1 + ln(2*pi)) + (1/2) * ln(sigmax^2)
% Entropy absnormal 2D = (2/2) * (1 + ln(2*pi)) + (1/2) * ln(det(covxy)) 
%-----------------------------------------------------------------------------------
% Entropy absnormal 2D = (2/2) * ln((2*pi)*exp(1)*(sigmax^2 * sigmay^2))^(1/2) 
% Entropy absnormal ND = (n/2) * ln((2*pi)*exp(1)*(sigmax^2 * sigmay^2 * ... * sigman^2))^(1/n) 
%-----------------------------------------------------------------------------------
%...................................................................................
ShannonEntropyTheoreticalCont002 = (1/2) * log((2*pi)*(logESDstd.^2)*exp(1)); %Okay.
ShannonEntropyTheoreticalCont003 = (1/2) * log((2*pi)*(logESDstd.^2)) + (1/2); %Okay.
ShannonEntropyTheoreticalCont004 = (1/2) * (1 + log(2*pi)) + (1/2)*log(logESDstd.^2); %Okay.
%...................................................................................
ShannonEntropyTheoreticalCont005 = log((2*pi)*(logESDstd.^2).*exp(logESDave + 1/2)); %WRONG!!!!
%...................................................................................
figure(11)
subplot(2,2,1)
imagesc(ShannonEntropyTheoreticalCont001)
title('ShannonEntropyTheoreticalCont -- 001')
colorbar_funhan(verticales)
subplot(2,2,2)
imagesc(ShannonEntropyTheoreticalCont002)
title('ShannonEntropyTheoreticalCont -- 002')
colorbar_funhan(verticales)
subplot(2,2,3)
imagesc(ShannonEntropyTheoreticalCont003)
title('ShannonEntropyTheoreticalCont -- 003')
colorbar_funhan(verticales)
subplot(2,2,4)
imagesc(ShannonEntropyTheoreticalCont004)
title('ShannonEntropyTheoreticalCont -- 004')
colorbar_funhan(verticales)
%...................................................................................
%===================================================================================
%...................................................................................
%%ShannonEntropyTheoreticalCont = ShannonEntropyTheoreticalCont001; %WEIRD...
ShannonEntropyTheoreticalCont = ShannonEntropyTheoreticalCont003; %THE GENERAL PATTERNS ARE OKAY BUT GIVES NEGATIVE VALUES.
%...................................................................................
%===================================================================================
%...................................................................................
%UXhdp = UXout(:,Jstep); %[d-1] 
%GXhdp = GXout(:,Jstep); %[d-1]
%...................................................................................
%UXssp = UXode(:,Jdays); %[d-1] 
%GXssp = GXode(:,Jdays); %[d-1] 
% Case with 2 traits (Le Gland, 17/07/2019)
% UXYhdp = UXYout(:,Jstep);
% GXYhdp = GXYout(:,Jstep);
%...................................................................................
UXYssp = UXYode(:,Jdays); 
GXYssp = GXYode(:,Jdays); 
%...................................................................................
SDINssp = SDINode(:,Jdays);
%...................................................................................
%FPHYTsspcontBis = UXssp .* PHYTsspcont;
FPHYTsspcontBis = UXYssp .* PHYTsspcont; % 2 traits (Le Gland, 17/07/2019)
%...................................................................................
%===================================================================================
% %...................................................................................
% FPHYThdpcont = FPHYToutcont(:,Jstep);
% EPHYThdpcont = EPHYToutcont(:,Jstep); 
% GPHYThdpcont = GPHYToutcont(:,Jstep); 
% MPHYThdpcont = MPHYToutcont(:,Jstep);
% %...................................................................................
% FZOOhdpcont = FZOOoutcont(:,Jstep);
% EZOOhdpcont = EZOOoutcont(:,Jstep);
% MZOOhdpcont = MZOOoutcont(:,Jstep);
% %...................................................................................
% FDINhdpcont = FDINoutcont(:,Jstep);
% FPONhdpcont = FPONoutcont(:,Jstep);
% %...................................................................................
%===================================================================================
% %...................................................................................
% FPHYTsspcont = FPHYTodecont(:,Jdays);
% EPHYTsspcont = EPHYTodecont(:,Jdays); 
% GPHYTsspcont = GPHYTodecont(:,Jdays); 
% MPHYTsspcont = MPHYTodecont(:,Jdays);
% %...................................................................................
% FZOOsspcont = FZOOodecont(:,Jdays);
% EZOOsspcont = EZOOodecont(:,Jdays);
% MZOOsspcont = MZOOodecont(:,Jdays);
% %...................................................................................
% FDINsspcont = FDINodecont(:,Jdays);
% FPONsspcont = FPONodecont(:,Jdays);
% %...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STATE VARIABLES DISCRETE TRAIT MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
[Vodedisc,todedisc] = jamstecrest_dailyAve(Voutdisc,deltat);
%...................................................................................
%===================================================================================
%...................................................................................
% % PHYhdpdisc3D = ones(ndepths,length(Jstep),nphy)*nan;
% PHYhdpdisc3D = ones(ndepths,length(Jstep),nxphy,nyphy)*nan; % 2 traits (Le Gland, 17/07/2019)
%...................................................................................
% PHYhdpdisc = Voutdisc(Jphy,Jstep); 
% ZOOhdpdisc = Voutdisc(Jzoo,Jstep);
% DINhdpdisc = Voutdisc(Jdin,Jstep);
% PONhdpdisc = Voutdisc(Jpon,Jstep);
%...................................................................................
%===================================================================================
%...................................................................................
%PHYodedisc3D = ones(ndepths,ndays*nyear,nphy)*nan;
PHYodedisc3D = ones(ndepths,ndays*nyear,nxphy,nyphy)*nan; % 2 traits (Le Gland, 17/07/2019)
%...................................................................................
PHYodedisc = Vodedisc(Jphy,:); 
ZOOodedisc = Vodedisc(Jzoo,:);
DINodedisc = Vodedisc(Jdin,:);
PONodedisc = Vodedisc(Jpon,:);
%...................................................................................
%===================================================================================
%...................................................................................
%PHYsspdisc3D = ones(ndepths,ndays,nphy)*nan;
PHYsspdisc3D = ones(ndepths,ndays,nxphy,nyphy)*nan; % 2 traits (Le Gland, 17/07/2019)
%...................................................................................
PHYsspdisc = PHYodedisc(:,Jdays);
ZOOsspdisc = ZOOodedisc(:,Jdays);
DINsspdisc = DINodedisc(:,Jdays);
PONsspdisc = PONodedisc(:,Jdays);
%...................................................................................
%===================================================================================
%...................................................................................
%uphysspdisc = uphyodedisc(:,Jdays,:); %size(depth,time,species) [d-1]
%gphysspdisc = gphyodedisc(:,Jdays,:); 
%...................................................................................
%FPHYsspdisc = FPHYodedisc(:,Jdays,:); 
%GPHYsspdisc = GPHYodedisc(:,Jdays,:); 
% Case with 2 traits (Le Gland, 17/07/2019)
% Comment to run without the discrete model (17/07/2020)
uphysspdisc = uphyodedisc(:,Jdays,:,:);
gphysspdisc = gphyodedisc(:,Jdays,:,:);
%uphysspdisc = zeros(ndepths,ndays,nxphy,nyphy);
%gphysspdisc = zeros(ndepths,ndays,nxphy,nyphy);
%...................................................................................
% FPHYsspdisc = FPHYodedisc(:,Jdays,:,:); 
% GPHYsspdisc = GPHYodedisc(:,Jdays,:,:); 
%...................................................................................
% $$$ uphyavessp = mean(uphysspdisc,3); %size(depth,time) [d-1] %DONT USE!!!
% $$$ gphyavessp = mean(gphysspdisc,3);
%...................................................................................
%===================================================================================
% %...................................................................................
% FPHYThdpdisc = FPHYToutdisc(:,Jstep);
% EPHYThdpdisc = EPHYToutdisc(:,Jstep); 
% GPHYThdpdisc = GPHYToutdisc(:,Jstep); 
% MPHYThdpdisc = MPHYToutdisc(:,Jstep);
% %...................................................................................
% FZOOhdpdisc = FZOOoutdisc(:,Jstep);
% EZOOhdpdisc = EZOOoutdisc(:,Jstep);
% MZOOhdpdisc = MZOOoutdisc(:,Jstep);
% %...................................................................................
% FDINhdpdisc = FDINoutdisc(:,Jstep);
% FPONhdpdisc = FPONoutdisc(:,Jstep);
% %...................................................................................
% %===================================================================================
% %...................................................................................
% FPHYTsspdisc = FPHYTodedisc(:,Jdays); 
% EPHYTsspdisc = EPHYTodedisc(:,Jdays); 
% GPHYTsspdisc = GPHYTodedisc(:,Jdays); 
% MPHYTsspdisc = MPHYTodedisc(:,Jdays);
% %...................................................................................
% FZOOsspdisc = FZOOodedisc(:,Jdays);
% EZOOsspdisc = EZOOodedisc(:,Jdays);
% MZOOsspdisc = MZOOodedisc(:,Jdays);
% %...................................................................................
% FDINsspdisc = FDINodedisc(:,Jdays);
% FPONsspdisc = FPONodedisc(:,Jdays);
% %...................................................................................
%===================================================================================
% $$$ %...................................................................................
% $$$ Izeros = find(PHYodedisc < sqrt(eps));
% $$$ Jzeros = find(PHYsspdisc < sqrt(eps));
% $$$ %...................................................................................
% $$$ PHYodedisc(Izeros) = 0d0;
% $$$ PHYsspdisc(Jzeros) = 0d0;
% $$$ %...................................................................................
%===================================================================================
% for iphy = 1:nphy
%     %....................................................................
%     Jphyi = Jphy(ndepths*(iphy-1)+1:ndepths*(iphy));
%     %....................................................................
%     iPHYhdp = PHYhdpdisc(Jphyi,:);
%     iPHYode = PHYodedisc(Jphyi,:);
%     iPHYssp = PHYsspdisc(Jphyi,:);
%     %....................................................................
%     PHYhdpdisc3D(:,:,iphy) = iPHYhdp; %Phyplankton [mmolP*m-3] (nphy,depths)
%     PHYodedisc3D(:,:,iphy) = iPHYode; %Phyplankton [mmolP*m-3] (nphy,depths)
%     PHYsspdisc3D(:,:,iphy) = iPHYssp; %Phyplankton [mmolP*m-3] (nphy,depths)
%     %....................................................................
%     numstrPHYhdp = ['PHY',sprintf('%02.0f',iphy),'out'];
%     numstrPHYode = ['PHY',sprintf('%02.0f',iphy),'ode'];
%     numstrPHYssp = ['PHY',sprintf('%02.0f',iphy),'ssp'];
%     %....................................................................
% % $$$     assignin('base',numstrPHYhdp,iPHYhdp) %NO USAR ESTA FORMA.
% % $$$     assignin('base',numstrPHYode,iPHYode)
% % $$$     assignin('base',numstrPHYssp,iPHYssp)
%     %....................................................................
%     myassign001 = [numstrPHYhdp,' = iPHYhdp;']; %USAR ESTA FORMA MEJOR FOR GLOBAL VARIABLES.
%     myassign002 = [numstrPHYode,' = iPHYode;']; %USAR ESTA FORMA MEJOR FOR GLOBAL VARIABLES.
%     myassign003 = [numstrPHYssp,' = iPHYssp;']; %USAR ESTA FORMA MEJOR FOR GLOBAL VARIABLES.
%     %....................................................................
%     eval(myassign001)
%     eval(myassign002)
%     eval(myassign003)
%     %....................................................................
% end
% Case with 2 traits (Le Gland, 17/07/2019)
for iphy = 1:nxphy
    for jphy = 1:nyphy 
        %....................................................................
        % WRONG ! Can cause errors when nxphy != nyphy
        % Jphyi = Jphy(ndepths*((jphy-1)*nyphy+(iphy-1))+1:ndepths*((jphy-1)*nyphy+iphy));
        Jphyi = Jphy(ndepths*((jphy-1)*nxphy+(iphy-1))+1:ndepths*((jphy-1)*nxphy+iphy));
        %....................................................................
        % iPHYhdp = PHYhdpdisc(Jphyi,:);
        %iPHYode = PHYodedisc(Jphyi,:);
        iPHYssp = PHYsspdisc(Jphyi,:);
        %....................................................................
        % PHYhdpdisc3D(:,:,iphy,jphy) = iPHYhdp; %Phyplankton [mmolP*m-3] (nphy,depths)
        %PHYodedisc3D(:,:,iphy,jphy) = iPHYode; %Phyplankton [mmolP*m-3] (nphy,depths)
        PHYsspdisc3D(:,:,iphy,jphy) = iPHYssp; %Phyplankton [mmolP*m-3] (nphy,depths)
        %....................................................................
        %numstrPHYhdp = ['PHY',sprintf('%02.0f',iphy),'out'];
        %numstrPHYode = ['PHY',sprintf('%02.0f',iphy),'ode'];
        numstrPHYssp = ['PHY',sprintf('%02.0f',iphy),'ssp'];
        %....................................................................
        %myassign001 = [numstrPHYhdp,' = iPHYhdp;']; %USAR ESTA FORMA MEJOR FOR GLOBAL VARIABLES.
        %myassign002 = [numstrPHYode,' = iPHYode;']; %USAR ESTA FORMA MEJOR FOR GLOBAL VARIABLES.
        myassign003 = [numstrPHYssp,' = iPHYssp;']; %USAR ESTA FORMA MEJOR FOR GLOBAL VARIABLES.
        %....................................................................
        %eval(myassign001)
        %eval(myassign002)
        eval(myassign003)
        %....................................................................
    end  
end
%%return
%........................................................................
% PHYThdpdisc = sum(PHYhdpdisc3D,3);
% PHYTodedisc = sum(PHYodedisc3D,3);
% PHYTsspdisc = sum(PHYsspdisc3D,3);
% Case with 2 traits (Le Gland, 17/07/2019)
% PHYThdpdisc = sum(PHYhdpdisc3D(:,:,:),3);
%PHYTodedisc = sum(PHYodedisc3D(:,:,:),3);
PHYTsspdisc = sum(PHYsspdisc3D(:,:,:),3);
%........................................................................
FPHYsspdisc3D = uphysspdisc .* PHYsspdisc3D; %[mmolN*m-3*d-1] 
GPHYsspdisc3D = gphysspdisc .* PHYsspdisc3D; %[mmolN*m-3*d-1] 
%...................................................................................
% FPHYTsspdisc = sum(FPHYsspdisc3D,3); %[mmolN*m-3*d-1] 
% GPHYTsspdisc = sum(GPHYsspdisc3D,3); %[mmolN*m-3*d-1] 
% Case with 2 traits (Le Gland, 17/07/2019)
FPHYTsspdisc = sum(FPHYsspdisc3D(:,:,:),3); 
GPHYTsspdisc = sum(GPHYsspdisc3D(:,:,:),3);
%...................................................................................
% uphytothdp = FPHYThdpdisc ./ PHYThdpdisc; %[d-1] 
% gphytothdp = GPHYThdpdisc ./ PHYThdpdisc; %[d-1] 
%...................................................................................
uphytotssp = FPHYTsspdisc ./ PHYTsspdisc; %[d-1] 
gphytotssp = GPHYTsspdisc ./ PHYTsspdisc; %[d-1] 
%...................................................................................
NTOTsspdisc = PHYTsspdisc + ZOOsspdisc + DINsspdisc + PONsspdisc;
%...................................................................................
%===================================================================================
%................................................................................... 
%logESDphy = xtraitrng; % Le Gland (05/06/2019), to be improved
logESDphy = xsizerng;
ESDphy = exp(logESDphy); %
[PHYsspcont3D] = jamstecrest_gaussecomodel1D_convertodiscrete(PHYTsspcont,logESDphy,logESDphysspAveCont,logESDphysspStdCont);
%................................................................................... 
% $$$ figure(5)
% $$$ plot(PHYsspcont3D(:),PHYsspdisc3D(:),'*')
% $$$ grid on
%................................................................................... 
%===================================================================================
%................................................................................... 
[logESDphyAve,logESDphyStd,B1coeff,B2coeff,R2coeff,ShannonEntropyDisc,ShannonEntropyTheoreticalDisc] = jamstecrest_discretemodel1D_lognormalcurvefit(PHYsspdisc3D,ESDphy,logESDphy,xsizedel);
%[logESDphyAve,logESDphyStd,B1coeff,B2coeff,R2coeff,ShannonEntropyDisc,ShannonEntropyTheoreticalDisc] = jamstecrest_discretemodel1D_lognormalcurvefit(PHYsspdisc3D,ESDphy,logESDphy,xtraitdel);
%...................................................................................
% This entropy is only about logESD and does not take temperature into account (Le Gland, 30/08/2019)
figure(112)
subplot(2,2,1)
imagesc(ShannonEntropyDisc)
title('ShannonEntropyDisc')
colorbar_funhan(verticales)
subplot(2,2,2)
imagesc(ShannonEntropyTheoreticalDisc)
title('ShannonEntropyTheoreticalDisc')
colorbar_funhan(verticales)
%...................................................................................
figure(113)
subplot(2,2,1)
plot(ShannonEntropyDisc(:),ShannonEntropyTheoreticalDisc(:),'*')
% $$$ set(gca,'Xlim',[0 0.5])
% $$$ set(gca,'Ylim',[-1 +1])
xlabel('ShannonEntropyDisc')
ylabel('ShannonEntropyTheoreticalDisc')
grid on
%...................................................................................
%===================================================================================
% $$$ %***********************************************************************************
% $$$ return

%%%%%%%%%%
%STOCKAGE:
%%%%%%%%%%
%===================================================================================
%...................................................................................
% ESDphyhdp3D = ones(ndepths,length(Jstep),nphy)*nan;
% ESDphyssp3D = ones(ndepths,ndays,nphy)*nan;
% Case with 2 traits (Le Gland, 18/07/2019)
% ESDphyhdp3D = ones(ndepths,length(Jstep),nxphy,nyphy)*nan;
% ESDphyssp3D = ones(ndepths,ndays,nxphy,nyphy)*nan;
% ESDphyhdp3D = ones(ndepths,length(Jstep),nxphy)*nan;
ESDphyssp3D = ones(ndepths,ndays,nxphy)*nan;
% SST case (Le Gland, 23/07/2019)
% ESDphyhdp3D = ones(ndepths,length(Jstep),nyphy)*nan;
% ESDphyssp3D = ones(ndepths,ndays,nyphy)*nan;

%...................................................................................
% for iphy = 1:nphy
%     esdphyi = ESDphy(iphy); 
%     ESDphyhdpi = esdphyi*ones(ndepths,length(Jstep));
%     ESDphysspi = esdphyi*ones(ndepths,ndays);
%     ESDphyhdp3D(:,:,iphy) = ESDphyhdpi; %[depth,time,species]
%     ESDphyssp3D(:,:,iphy) = ESDphysspi; %[depth,time,species]
% end
% Case with 2 traits (Le Gland, 18/07/2019)
for iphy = 1:nxphy
%for iphy = 1:nyphy % SST case (Le Gland, 23/07/2019)
    esdphyi = ESDphy(iphy); 
    % ESDphyhdpi = esdphyi*ones(ndepths,length(Jstep));
    ESDphysspi = esdphyi*ones(ndepths,ndays);
    % ESDphyhdp3D(:,:,iphy) = ESDphyhdpi; %[depth,time,species]
    ESDphyssp3D(:,:,iphy) = ESDphysspi; %[depth,time,species]
%     for jphy = 1:nyphy
%         % 4D arrays which values depend only on iphy
%         ESDphyhdp3D(:,:,iphy,jphy) = ESDphyhdpi; %[depth,time,species]
%         ESDphyssp3D(:,:,iphy,jphy) = ESDphysspi; %[depth,time,species]
%     end
end
%...................................................................................
%===================================================================================
%...................................................................................
% logESDphyhdp3D = ones(ndepths,length(Jstep),nphy)*nan;
% logESDphyssp3D = ones(ndepths,ndays,nphy)*nan;
% Case with 2 traits (Le Gland, 18/07/2019)
% logESDphyhdp3D = ones(ndepths,length(Jstep),nxphy,nyphy)*nan;
% logESDphyssp3D = ones(ndepths,ndays,nxphy,nyphy)*nan;
% logESDphyhdp3D = ones(ndepths,length(Jstep),nxphy)*nan;
logESDphyssp3D = ones(ndepths,ndays,nxphy)*nan;
% logESDphyhdp3D = ones(ndepths,length(Jstep),nyphy)*nan;
% logESDphyssp3D = ones(ndepths,ndays,nyphy)*nan;
%...................................................................................
% for iphy = 1:nphy
%     logesdphyi = logESDphy(iphy); 
%     logESDphyhdpi = logesdphyi*ones(ndepths,length(Jstep));
%     logESDphysspi = logesdphyi*ones(ndepths,ndays);
%     logESDphyhdp3D(:,:,iphy) = logESDphyhdpi; %[depth,time,species]
%     logESDphyssp3D(:,:,iphy) = logESDphysspi; %[depth,time,species]
% end
% Case with 2 traits (Le Gland, 18/07/2019)
for iphy = 1:nxphy
    logesdphyi = logESDphy(iphy);
%for iphy = 1:nyphy 
%    logesdphyi = ytrait(iphy); % SST (Le Gland, 23/07/2019) 
    % logESDphyhdpi = logesdphyi*ones(ndepths,length(Jstep));
    logESDphysspi = logesdphyi*ones(ndepths,ndays);
    % logESDphyhdp3D(:,:,iphy) = logESDphyhdpi; %[depth,time,species]
    logESDphyssp3D(:,:,iphy) = logESDphysspi; %[depth,time,species]
%     for jphy = 1:nyphy
%         logESDphyhdp3D(:,:,iphy,jphy) = logESDphyhdpi; %[depth,time,species]
%         logESDphyssp3D(:,:,iphy,jphy) = logESDphysspi; %[depth,time,species]
%     end
end
%...................................................................................
%===================================================================================
%...................................................................................

% Intent to visualize both logESD and Topt (Le Gland, 02/09/2019)
% TOPTphyhdp3D = ones(ndepths,length(Jstep),nyphy)*nan;
TOPTphyssp3D = ones(ndepths,ndays,nyphy)*nan;
for iphy = 1:nyphy 
    Toptphyi = ytrait(iphy); % SST (Le Gland, 23/07/2019) 
    % TOPTphyhdpi = Toptphyi*ones(ndepths,length(Jstep));
    TOPTphysspi = Toptphyi*ones(ndepths,ndays);
    % TOPTphyhdp3D(:,:,iphy) = TOPTphyhdpi; %[depth,time,species]
    TOPTphyssp3D(:,:,iphy) = TOPTphysspi; %[depth,time,species]
end

%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%APPROXIMATE THE CONTINUOUS DISTRIBUTION STATISTICS FROM THE DISCRETE SIMULATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%[ESDphyhdpAveDiscBis,ESDphyhdpStdDiscBis] = jamstecrest_geometricmean(PHYhdpdisc3D,ESDphyhdp3D); 
%[ESDphysspAveDiscBis,ESDphysspStdDiscBis] = jamstecrest_geometricmean(PHYsspdisc3D,ESDphyssp3D); 
%...................................................................................
%[logESDphyhdpAveDisc,logESDphyhdpStdDisc] = jamstecrest_geometricmean(PHYhdpdisc3D,logESDphyhdp3D); 
%[logESDphysspAveDisc,logESDphysspStdDisc] = jamstecrest_geometricmean(PHYsspdisc3D,logESDphyssp3D); 
%...................................................................................
%[TOPTphyhdpAveDisc,TOPTphyhdpStdDisc] = jamstecrest_geometricmean(PHYhdpdisc3D,TOPTphyhdp3D); 
%[TOPTphysspAveDisc,TOPTphysspStdDisc] = jamstecrest_geometricmean(PHYsspdisc3D,TOPTphyssp3D);
%...................................................................................
% Correlation of the discrete model (Le Gland, 03/09/2019)
% By the way, I suspect geometricmean contains errors and could be cleaned
%[logESDphyhdpAveDisc,logESDphyhdpStdDisc,TOPTphyhdpAveDisc,TOPTphyhdpStdDisc,phyhdpCorDisc] = jamstecrest_covariance(PHYhdpdisc3D,logESDphyhdp3D,TOPTphyhdp3D);
[logESDphysspAveDisc,logESDphysspStdDisc,TOPTphysspAveDisc,TOPTphysspStdDisc,physspCorDisc] = jamstecrest_covariance(PHYsspdisc3D,logESDphyssp3D,TOPTphyssp3D);
%ESDphyhdpAveDiscBis = exp(logESDphyhdpStdDisc);
ESDphysspAveDiscBis = exp(logESDphysspStdDisc);
% Entropy of the discrete model (Le Gland, 05/11/2019)
%phyhdpEntDisc = 2*pi*exp(1).*logESDphyhdpStdDisc.*TOPTphyhdpStdDisc.*sqrt(1-phyhdpCorDisc.^2);
physspEntDisc = 2*pi*exp(1).*logESDphysspStdDisc.*TOPTphysspStdDisc.*sqrt(1-physspCorDisc.^2);
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SHANNON CONTINOUS ENTROPY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
% YOPT is a confusing name, since is applies to logESD
% I do not change it for now, since this part is not important (Le Gland, 02/09/2019)
%................................................................................... 
XOPTave = logESDphysspAveCont + 1;
XOPTstd = logESDphysspStdCont + 1;
%................................................................................... 
YOPTave = logESDphysspAveDisc + 1;
YOPTstd = logESDphysspStdDisc + 1;
%................................................................................... 
%===================================================================================
%NORMALIZING VALUES BETWEEN 0 AND 1:
%................................................................................... 
XOPTave_Star = XOPTave / max(XOPTave(:));
XOPTstd_Star = XOPTstd / max(XOPTstd(:));
%................................................................................... 
YOPTave_Star = YOPTave / max(YOPTave(:));
YOPTstd_Star = YOPTstd / max(YOPTstd(:));
%................................................................................... 
%===================================================================================
%................................................................................... 

fignum = 101;
%................................................................................... 
[ShannonEntropyTheoreticalAbsnormal_Cont,ShannonEntropyTheoreticalLognormal_Cont] = jamstecrest_ShannonContinuousEntropy1D(XOPTave_Star,XOPTstd_Star,fignum,mypackages);
%................................................................................... 
fignum = 102;
%................................................................................... 
[ShannonEntropyTheoreticalAbsnormal_Disc,ShannonEntropyTheoreticalLognormal_Disc] = jamstecrest_ShannonContinuousEntropy1D(YOPTave_Star,YOPTstd_Star,fignum,mypackages);
%...................................................................................
Amin_Absnormal = min([ShannonEntropyTheoreticalAbsnormal_Cont(:);ShannonEntropyTheoreticalAbsnormal_Disc(:)]);
Amin_Lognormal = min([ShannonEntropyTheoreticalLognormal_Cont(:);ShannonEntropyTheoreticalLognormal_Disc(:)]);
%................................................................................... 
Amax_Absnormal = max([ShannonEntropyTheoreticalAbsnormal_Cont(:);ShannonEntropyTheoreticalAbsnormal_Disc(:)]);
Amax_Lognormal = max([ShannonEntropyTheoreticalLognormal_Cont(:);ShannonEntropyTheoreticalLognormal_Disc(:)]);
%................................................................................... 
figure(105)
subplot(2,2,1)
imagesc(ShannonEntropyTheoreticalAbsnormal_Cont,[Amin_Absnormal Amax_Absnormal])
title('Entropy of Abs Normal Distribution -- Cont')
colorbar_funhan(verticales)
subplot(2,2,2)
imagesc(ShannonEntropyTheoreticalLognormal_Cont,[Amin_Lognormal Amax_Lognormal])
title('Entropy of Log Normal Distribution -- Cont')
colorbar_funhan(verticales)
grid on
subplot(2,2,3)
imagesc(ShannonEntropyTheoreticalAbsnormal_Disc,[Amin_Absnormal Amax_Absnormal])
title('Entropy of Abs Normal Distribution -- Disc')
colorbar_funhan(verticales)
grid on
subplot(2,2,4)
imagesc(ShannonEntropyTheoreticalLognormal_Disc,[Amin_Lognormal Amax_Lognormal])
title('Entropy of Log Normal Distribution -- Disc')
colorbar_funhan(verticales)
grid on
%................................................................................... 
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE UPPER AND LOWER VALUES OF THE CONTINUOUS DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%COMPUTING AVERAGE AND RANGE OF CELL SIZE IN ABSOLUTE SCALE FROM CELL SIZE:
%(NOTE: I DON'T THINK THIS IS THE RIGHT WAY TO DO IT - DONT USE) 
%...................................................................................
%ESDphysspAveUpperDiscBis = ESDphysspAveDiscBis + ESDphysspStdDiscBis; %WRONG!!!
%ESDphysspAveLowerDiscBis = ESDphysspAveDiscBis - ESDphysspStdDiscBis; 
%...................................................................................
%ESDphysspCVdiscBis = ESDphysspStdDiscBis ./ ESDphysspAveDiscBis; %WRONG!!!
%...................................................................................
%===================================================================================
%...................................................................................
%[ESDphyhdpAveUpperCont,ESDphyhdpAveLowerCont,ESDphyhdpAveCont,ESDphyhdpStdCont,ESDphyhdpCVcont] = jamstecrest_meansigmaESD(logESDphyhdpAveCont,logESDphyhdpStdCont);
[ESDphysspAveUpperCont,ESDphysspAveLowerCont,ESDphysspAveCont,ESDphysspStdCont,ESDphysspCVcont] = jamstecrest_meansigmaESD(logESDphysspAveCont,logESDphysspStdCont);
%...................................................................................
%[ESDphyhdpAveUpperDisc,ESDphyhdpAveLowerDisc,ESDphyhdpAveDisc,ESDphyhdpStdDisc,ESDphyhdpCVdisc] = jamstecrest_meansigmaESD(logESDphyhdpAveDisc,logESDphyhdpStdDisc);
[ESDphysspAveUpperDisc,ESDphysspAveLowerDisc,ESDphysspAveDisc,ESDphysspStdDisc,ESDphysspCVdisc] = jamstecrest_meansigmaESD(logESDphysspAveDisc,logESDphysspStdDisc);
%...................................................................................
%[ESDphyhdpAveCont002,ESDphyhdpVarCont002] = funlognstat(logESDphyhdpAveCont,logESDphyhdpStdCont);
[ESDphysspAveCont002,ESDphysspVarCont002] = funlognstat(logESDphysspAveCont,logESDphysspStdCont);
%...................................................................................
%[ESDphyhdpAveDisc002,ESDphyhdpVarDisc002] = funlognstat(logESDphyhdpAveDisc,logESDphyhdpStdDisc);
[ESDphysspAveDisc002,ESDphysspVarDisc002] = funlognstat(logESDphysspAveDisc,logESDphysspStdDisc);
%...................................................................................
% $$$ [ESDphyhdpAveCont002,ESDphyhdpVarCont002] = funlognstat(logESDphyhdpAveCont,sqrt(logESDphyhdpStdCont));
% $$$ [ESDphysspAveCont002,ESDphysspVarCont002] = funlognstat(logESDphysspAveCont,sqrt(logESDphysspStdCont));
% $$$ %...................................................................................
% $$$ [ESDphyhdpAveDisc002,ESDphyhdpVarDisc002] = funlognstat(logESDphyhdpAveDisc,sqrt(logESDphyhdpStdDisc));
% $$$ [ESDphysspAveDisc002,ESDphysspVarDisc002] = funlognstat(logESDphysspAveDisc,sqrt(logESDphysspStdDisc));
%...................................................................................
%ESDphyhdpStdCont002 = sqrt(ESDphyhdpVarCont002);
ESDphysspStdCont002 = sqrt(ESDphysspVarCont002);
%...................................................................................
%ESDphyhdpStdDisc002 = sqrt(ESDphyhdpVarDisc002);
ESDphysspStdDisc002 = sqrt(ESDphysspVarDisc002);
%...................................................................................
%===================================================================================
%...................................................................................
ESDphysspDelta_Cont = (ESDphysspAveUpperCont - ESDphysspAveLowerCont);
ESDphysspDelta_Disc = (ESDphysspAveUpperDisc - ESDphysspAveLowerDisc);
%...................................................................................
%===================================================================================
%...................................................................................
%figure(400)
%subplot(2,2,1)
%plot(ESDphyhdpAveCont(:),ESDphyhdpAveCont002(:),'*')
%title('ESD ave - Cont -- HPD')
%subplot(2,2,2)
%plot(ESDphysspAveCont(:),ESDphysspAveCont002(:),'*')
%title('ESD ave - Cont -- SSP')
%subplot(2,2,3)
%plot(ESDphyhdpAveDisc(:),ESDphyhdpAveDisc002(:),'*')
%title('ESD ave - Disc -- HPD')
%subplot(2,2,4)
%plot(ESDphysspAveDisc(:),ESDphysspAveDisc002(:),'*')
%title('ESD ave - Disc -- SSP')
%...................................................................................
%figure(410)
%subplot(2,2,1)
%plot(ESDphyhdpStdCont(:),ESDphyhdpStdCont002(:),'*')
%title('ESD std - Cont -- HPD')
%subplot(2,2,2)
%plot(ESDphysspStdCont(:),ESDphysspStdCont002(:),'*')
%title('ESD std - Cont -- SSP')
%subplot(2,2,3)
%plot(ESDphyhdpStdDisc(:),ESDphyhdpStdDisc002(:),'*')
%title('ESD std - Disc -- HPD')
%subplot(2,2,4)
%plot(ESDphysspStdDisc(:),ESDphysspStdDisc002(:),'*')
%title('ESD std - Disc -- SSP')
%%return
%...................................................................................
%figure(500)
%%subplot(2,4,1)
%subplot(5,4,1)
%imagesc(logESDphyhdpAveCont)
%title('logESD ave - cont (hdp)')
%colorbar_funhan(horizontal)
%%subplot(2,4,2)
%subplot(5,4,2)
%imagesc(logESDphysspAveCont)
%title('logESD ave - cont (ssp)')
%colorbar_funhan(horizontal)
%%subplot(2,4,3)
%subplot(5,4,3)
%imagesc(logESDphyhdpStdCont)
%title('logESD std - cont (hdp)')
%colorbar_funhan(horizontal)
%%subplot(2,4,4)
%subplot(5,4,4)
%imagesc(logESDphysspStdCont)
%title('logESD std - cont (ssp)')
%colorbar_funhan(horizontal)
%%subplot(2,4,1+4)
%subplot(5,4,1+4)
%imagesc(logESDphyhdpAveDisc)
%title('logESD ave - disc (hdp)')
%colorbar_funhan(horizontal)
%%subplot(2,4,2+4)
%subplot(5,4,2+4)
%imagesc(logESDphysspAveDisc)
%title('logESD ave - disc (ssp)')
%colorbar_funhan(horizontal)
%%subplot(2,4,3+4)
%subplot(5,4,3+4)
%imagesc(logESDphyhdpStdDisc)
%title('logESD std - disc (hdp)')
%colorbar_funhan(horizontal)
%%subplot(2,4,4+4)
%subplot(5,4,4+4)
%imagesc(logESDphysspStdDisc)
%title('logESD std - disc (ssp)')
%colorbar_funhan(horizontal)

%% SST part (Le Gland, 02/09/2019)
%subplot(5,4,1+8)
%imagesc(TOPTphyhdpAveCont)
%title('Topt ave - cont (hdp)')
%colorbar_funhan(horizontal)
%subplot(5,4,2+8)
%imagesc(TOPTphysspAveCont)
%title('Topt ave - cont (ssp)')
%colorbar_funhan(horizontal)
%subplot(5,4,3+8)
%imagesc(TOPTphyhdpStdCont)
%title('Topt std - cont (hdp)')
%colorbar_funhan(horizontal)
%subplot(5,4,4+8)
%imagesc(TOPTphysspStdCont)
%title('Topt std - cont (ssp)')
%colorbar_funhan(horizontal)
%subplot(5,4,1+12)
%imagesc(TOPTphyhdpAveDisc)
%title('Topt ave - disc (hdp)')
%colorbar_funhan(horizontal)
%subplot(5,4,2+12)
%imagesc(TOPTphysspAveDisc)
%title('Topt ave - disc (ssp)')
%colorbar_funhan(horizontal)
%subplot(5,4,3+12)
%imagesc(TOPTphyhdpStdDisc)
%title('Topt std - disc (hdp)')
%colorbar_funhan(horizontal)
%subplot(5,4,4+12)
%imagesc(TOPTphysspStdDisc)
%title('Topt std - disc (ssp)')
%colorbar_funhan(horizontal)

%%Correlation part (Le Gland, 02/09/2019)
%subplot(5,4,1+16)
%imagesc(phyhdpCorCont)
%title('Topt/logESD cor - cont (hdp)')
%colorbar_funhan(horizontal)
%subplot(5,4,2+16)
%imagesc(physspCorCont)
%title('Topt/logESD cor - cont (ssp)')
%colorbar_funhan(horizontal)
%subplot(5,4,3+16)
%imagesc(phyhdpCorDisc)
%title('Topt/logESD cor - disc (hdp)')
%colorbar_funhan(horizontal)
%subplot(5,4,4+16)
%imagesc(physspCorDisc)
%title('Topt/logESD cor - disc (ssp)')
%colorbar_funhan(horizontal)

% Print this figure too (Le Gland, 17/09/2019)
%%print('-dpng','-r300','jamstecrest_gaussecomodel1D_fig500.png')

%...................................................................................
%figure(600)
%subplot(2,4,1)
%imagesc(ESDphyhdpAveCont)
%title('ESD ave - cont (hdp)')
%colorbar_funhan(horizontal)
%subplot(2,4,2)
%imagesc(ESDphysspAveCont)
%title('ESD ave - cont (ssp)')
%colorbar_funhan(horizontal)
%subplot(2,4,3)
%imagesc(ESDphyhdpStdCont)
%title('ESD std - cont (hdp)')
%colorbar_funhan(horizontal)
%subplot(2,4,4)
%imagesc(ESDphysspStdCont)
%title('ESD std - cont (ssp)')
%colorbar_funhan(horizontal)
%subplot(2,4,1+4)
%imagesc(ESDphyhdpAveDisc)
%title('ESD ave - disc (hdp)')
%colorbar_funhan(horizontal)
%subplot(2,4,2+4)
%imagesc(ESDphysspAveDisc)
%title('ESD ave - disc (ssp)')
%colorbar_funhan(horizontal)
%subplot(2,4,3+4)
%imagesc(ESDphyhdpStdDisc)
%title('ESD std - disc (hdp)')
%colorbar_funhan(horizontal)
%subplot(2,4,4+4)
%imagesc(ESDphysspStdDisc)
%title('ESD std - disc (ssp)')
%colorbar_funhan(horizontal)
%...................................................................................
%figure(700)
%subplot(2,4,1)
%imagesc(ESDphyhdpAveCont002)
%title('ESD ave Bis - cont (hdp)')
%colorbar_funhan(horizontal)
%subplot(2,4,2)
%imagesc(ESDphysspAveCont002)
%title('ESD ave Bis - cont (ssp)')
%colorbar_funhan(horizontal)
%subplot(2,4,3)
%imagesc(ESDphyhdpStdCont002)
%title('ESD std Bis - cont (hdp)')
%colorbar_funhan(horizontal)
%subplot(2,4,4)
%imagesc(ESDphysspStdCont002)
%title('ESD std Bis - cont (ssp)')
%colorbar_funhan(horizontal)
%subplot(2,4,1+4)
%imagesc(ESDphyhdpAveDisc002)
%title('ESD ave Bis - disc (hdp)')
%colorbar_funhan(horizontal)
%subplot(2,4,2+4)
%imagesc(ESDphysspAveDisc002)
%title('ESD ave Bis - disc (ssp)')
%colorbar_funhan(horizontal)
%subplot(2,4,3+4)
%imagesc(ESDphyhdpStdDisc002)
%title('ESD std Bis - disc (hdp)')
%colorbar_funhan(horizontal)
%subplot(2,4,4+4)
%imagesc(ESDphysspStdDisc002)
%title('ESD std Bis - disc (ssp)')
%colorbar_funhan(horizontal)
%...................................................................................
%===================================================================================
%...................................................................................
%figure(40)
%subplot(2,2,1) 
%plot(ESDphysspAveDisc(:),ESDphysspAveDiscBis(:),'*') 
%axis([0 5, 0 5])
%xlabel('from linear scale [um]')
%ylabel('from lognep scale [um]')
%title('ESD ave')
%grid on 
%subplot(2,2,2) 
%plot(ESDphysspStdDisc(:),ESDphysspStdDiscBis(:),'*') 
%axis([0 5, 0 5])
%xlabel('from linear scale [um]')
%ylabel('from lognep scale [um]')
%title('ESD std')
%grid on 
%subplot(2,2,3) 
%plot(ESDphysspCVdisc(:),ESDphysspCVdiscBis(:),'*')
%axis([0 1, 0 1])
%xlabel('from linear scale [n.d.]')
%ylabel('from lognep scale [n.d.]')
%title('ESD cv')
%grid on 
%...................................................................................
%===================================================================================
%figure(50) %Should be the same!!!
%subplot(2,2,1)
%imagesc(logESDphyAve)
%colorbar
%subplot(2,2,2)
%imagesc(logESDphyStd)
%colorbar
%subplot(2,2,3)
%imagesc(logESDphysspAveDisc) 
%colorbar
%subplot(2,2,4)
%imagesc(logESDphysspStdDisc) 
%colorbar
%...................................................................................
%figure(60) %Should be: B1 = 0.0 (intercept), B2 = 1.0 (slope), R2 = 1.0 (coeff.det)
%subplot(2,2,1)
%imagesc(B1coeff) 
%caxis([0.0 + [-0.5, +0.5]])
%colorbar
%subplot(2,2,2)
%imagesc(B2coeff) 
%caxis([1.0 + [-0.5, +0.5]])
%colorbar
%subplot(2,2,3)
%imagesc(R2coeff)
%caxis([0 1])
%colorbar
%...................................................................................
%pause(1)
%close(40)
%close(50)
%close(60)
%===================================================================================
%...................................................................................
[HindexCont,JindexCont,SindexCont,EindexCont,DindexCont] = jamstecrest_gaussecomodel1D_ShannonIndex(PHYsspcont3D,mypackages);
%...................................................................................
[HindexDisc,JindexDisc,SindexDisc,EindexDisc,DindexDisc] = jamstecrest_gaussecomodel1D_ShannonIndex(PHYsspdisc3D,mypackages);
%...................................................................................
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% $$$ %...................................................................................
% $$$ logbase = exp(1);
% $$$ for jdepth = 1:ndepths
% $$$     jphy = squeeze(PHYsspcont3D(jdepth,:,:));
% $$$     [hindex,junkarg] = myShannonWienerIndex(jphy',logbase);
% $$$     HindexContBis(jdepth,:) = hindex;
% $$$ end
% $$$ %...................................................................................
% $$$ figure(3)
% $$$ plot(HindexCont(:),HindexContBis(:),'*')
% $$$ axis([0 4, 0 4])
% $$$ grid on
% $$$ %...................................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%===================================================================================
%HIGH PASS FILTER FOR MASKING DIVERSITY ESTIMATES WHERE THERE IS ALMOST NO PHY BIOMASS:
%...................................................................................
%PHY = PHYsspdisc3D;
%PHY = sum(PHYsspdisc3D,4); % 2 traits (Le Gland, 18/07/2019) 
%PHY = sum(PHYsspdisc3D,3); % SST case (Le Gland, 23/07/2019)
%sumPHY = sum(PHY,3);

%PHYsspdisc3Dpcnt = PHY ./ repmat(sumPHY,[1,1,nphy]);
%PHYsspdisc3Dpcnt = PHY ./ repmat(sumPHY,[1,1,nxphy]); % 2 traits (Le Gland, 19/07/2019)

%Symmetrical treatment of the 2 traits (Le Gland, 02/09/2019)
PHY = PHYsspdisc3D;
sumPHY = sum(sum(PHY,4),3); 
PHYsspdisc3Dpcnt = PHY ./ repmat(sumPHY,[1,1,nxphy,nyphy]);

%...................................................................................
% $$$ kphy = median(sumPHY(:));
% $$$ kphy = 0.01*mean(sumPHY(:));
kphy = 0.05*max(sumPHY(:));
% $$$ kphy = 0.01*max(sumPHY(:));
% $$$ kphy = 0.001*max(sumPHY(:));
%...................................................................................
HighPassFilterMaskPhydisc = (sumPHY.^2) ./ (kphy^2 + sumPHY.^2);
%...................................................................................
Iones = find(HighPassFilterMaskPhydisc < 0.6);
%...................................................................................
HighPassFilterMaskPhydisc(Iones) = nan;
%...................................................................................
HindexMasked = HindexDisc .* HighPassFilterMaskPhydisc;
JindexMasked = JindexDisc .* HighPassFilterMaskPhydisc;
SindexMasked = SindexDisc .* HighPassFilterMaskPhydisc;
EindexMasked = EindexDisc .* HighPassFilterMaskPhydisc;
DindexMasked = DindexDisc .* HighPassFilterMaskPhydisc;
%...................................................................................
%===================================================================================
%...................................................................................
% $$$ figure(2400)
% $$$ subplot(2,2,1)
% $$$ semilogx(sumPHY(:),HighPassFilterMaskPhydisc(:),'*')
% $$$ xlabel('PHYT')
% $$$ ylabel('High Pass Filter')
% $$$ grid on
% $$$ subplot(2,2,2)
% $$$ imagesc(HighPassFilterMaskPhydisc,[0 1])
% $$$ colorbar_funhan(verticales)
% $$$ grid on
%...................................................................................
maxSindex = max([max(SindexCont(:)),max(SindexDisc(:))]);
minSindex = min([min(SindexCont(:)),min(SindexDisc(:))]);
%...................................................................................
maxEindex =  ceil(max([max(EindexCont(:)),max(EindexDisc(:))]));
minEindex = floor(min([min(EindexCont(:)),min(EindexDisc(:))]));
%...................................................................................
maxShannonEntropy = max([max(ShannonEntropyTheoreticalCont(:)),max(ShannonEntropyTheoreticalDisc(:))]);
minShannonEntropy = min([min(ShannonEntropyTheoreticalCont(:)),min(ShannonEntropyTheoreticalDisc(:))]);
%...................................................................................
figure(2500)
hplot = subplot(2,3,1);
imagesc(ShannonEntropyTheoreticalCont)
colorbar_funhan(verticales)
title(hplot,'Shannon Entropy - Cont')
grid on
hplot = subplot(2,3,2);
imagesc(SindexCont,[minSindex maxSindex])
colorbar_funhan(verticales)
title(hplot,'Species Richness - Cont')
grid on
hplot = subplot(2,3,3);
imagesc(EindexCont,[minEindex maxEindex])
colorbar_funhan(verticales)
title(hplot,'exp (Shannon index) - Cont')
grid on
hplot = subplot(2,3,1+3);
imagesc(ShannonEntropyTheoreticalDisc)
colorbar_funhan(verticales)
title(hplot,'Shannon Entropy - Disc')
grid on
hplot = subplot(2,3,2+3);
imagesc(SindexDisc,[minSindex maxSindex])
colorbar_funhan(verticales)
title(hplot,'Species Richness - Disc')
grid on
hplot = subplot(2,3,3+3);
imagesc(EindexDisc,[minEindex maxEindex])
colorbar_funhan(verticales)
title(hplot,'exp (Shannon index) - Disc')
grid on
%...................................................................................
%%print('-dpng','-r300','jamstecrest_discretemodel1D_fig003.png')
%...................................................................................
% $$$ figure(2005)
% $$$ hplot = subplot(2,2,1);
% $$$ plot(exp(ShannonEntropyDisc(:)),ESDphysspCVdisc(:),'*')
% $$$ grid on
%...................................................................................
% $$$ figure(2510)
% $$$ hplot = subplot(2,2,1);
% $$$ mynanimagesc(DindexMasked)
% $$$ title(hplot,'Simpson index')
% $$$ grid on
% $$$ hplot = subplot(2,2,2);
% $$$ mynanimagesc(JindexMasked)
% $$$ title(hplot,'Shannon index normalized [0 - 1]')
% $$$ grid on
% $$$ hplot = subplot(2,2,3);
% $$$ mynanimagesc(SindexMasked)
% $$$ title(hplot,'Species Richness')
% $$$ grid on
% $$$ hplot = subplot(2,2,4);
% $$$ mynanimagesc(EindexMasked)
% $$$ title(hplot,'exp (Shannon index)')
% $$$ grid on
%...................................................................................
% $$$ print('-dpng','-r300','jamstecrest_discretemodel1D_fig004.png')
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPECIFIC PRODUCTION RATES (PER DAY):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
% MUPhdpcont = FPHYThdpcont ./ PHYThdpcont; %[d-1]
% MUPodecont = FPHYTodecont ./ PHYTodecont; %[d-1] 
% MUPsspcont = FPHYTsspcont ./ PHYTsspcont; %[d-1]
% Use only Voutcont variables (Le Gland, 27/11/2019)
% MUPhdpcont = UXYhdp;
% MUPodecont = UXYode;
MUPsspcont = UXYssp;
%...................................................................................
% MUPhdpdisc = FPHYThdpdisc ./ PHYThdpdisc; %[d-1]
% MUPodedisc = FPHYTodedisc ./ PHYTodedisc; %[d-1] 
% MUPsspdisc = FPHYTsspdisc ./ PHYTsspdisc; %[d-1] 
% Use only Voutcont variables (Le Gland, 27/11/2019)
% MUPhdpdisc is unavailable: only daily outputs from discrete model
% MUPodedisc = sum(uphyodedisc(:,:,:).*PHYodedisc3D(:,:,:),3)./PHYTodedisc;
MUPsspdisc = sum(uphysspdisc(:,:,:).*PHYsspdisc3D(:,:,:),3)./PHYTsspdisc;
%MUPodedisc = 0;
%MUPsspdisc = 0;

%...................................................................................
% $$$ MUPodecontBis = UX; 
% $$$ MUPsspcontBis = UXssp; 
%...................................................................................
%===================================================================================
%...................................................................................
% MUZhdpcont = FZOOhdpcont ./ ZOOhdpcont;  %[d-1]
% MUZodecont = FZOOodecont ./ ZOOodecont;  %[d-1]
% MUZsspcont = FZOOsspcont ./ ZOOsspcont;  %[d-1]
% Use only Voutcont variables (Le Gland, 27/11/2019)
% MUZhdpcont = GXYhdp .* (PHYThdpcont ./ ZOOhdpcont);
% MUZodecont = GXYode .* (PHYTodecont ./ ZOOodecont);
MUZsspcont = GXYssp .* (PHYTsspcont ./ ZOOsspcont);
%...................................................................................
% MUZhdpdisc = FZOOhdpdisc ./ ZOOhdpdisc;  %[d-1]
% MUZodedisc = FZOOodedisc ./ ZOOodedisc;  %[d-1]
% MUZsspdisc = FZOOsspdisc ./ ZOOsspdisc;  %[d-1]
% Use only Voutcont variables (Le Gland, 27/11/2019)
% MUZhdpdisc is unavailable: only daily outputs from discrete model
% MUZodedisc = sum(gphyodedisc(:,:,:).*PHYodedisc3D(:,:,:),3)./ZOOodedisc;  %
MUZsspdisc = sum(gphysspdisc(:,:,:).*PHYsspdisc3D(:,:,:),3)./ZOOsspdisc;  %
% MUZodedisc = 0;
%MUZsspdisc = 0;
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%
%SAVE MATLAB BINARY:
%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
save('jamstecrest_gaussecomodel1D_workspace.mat') 
%...................................................................................
%===================================================================================
%***********************************************************************************
%%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTS FOR CONTINOUS AND DISCRETE MODELS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
minESDphyAveLower = min([min(ESDphysspAveLowerDisc(:)),min(ESDphysspAveLowerCont(:))])
maxESDphyAveLower = max([max(ESDphysspAveLowerDisc(:)),max(ESDphysspAveLowerCont(:))])
%...................................................................................
minESDphyAveUpper = min([min(ESDphysspAveUpperDisc(:)),min(ESDphysspAveUpperCont(:))])
maxESDphyAveUpper = max([max(ESDphysspAveUpperDisc(:)),max(ESDphysspAveUpperCont(:))])
%...................................................................................
minESDphyAve = min([min(ESDphysspAveDisc(:)),min(ESDphysspAveCont(:))])
maxESDphyAve = max([max(ESDphysspAveDisc(:)),max(ESDphysspAveCont(:))])
%...................................................................................
minESDphyAveLower = -eps;
minESDphyAveUpper = -eps;
minESDphyAve = -eps;
%...................................................................................
figure(3010)
hplot = subplot(2,3,1);
imagesc(ESDphysspAveLowerCont,[minESDphyAveLower maxESDphyAveLower])
colorbar_funhan(verticales)
title(hplot,'ESD lower - cont')
grid on
hplot = subplot(2,3,2);
imagesc(ESDphysspAveCont,[minESDphyAve maxESDphyAve])
colorbar_funhan(verticales)
title(hplot,'ESD ave - cont')
grid on
hplot = subplot(2,3,3);
imagesc(ESDphysspAveUpperCont,[minESDphyAveUpper maxESDphyAveUpper])
colorbar_funhan(verticales)
title(hplot,'ESD upper - cont')
grid on
%...................................................................................
hplot = subplot(2,3,1+3);
imagesc(ESDphysspAveLowerDisc,[minESDphyAveLower maxESDphyAveLower])
colorbar_funhan(verticales)
title(hplot,'ESD lower - disc')
grid on
hplot = subplot(2,3,2+3);
imagesc(ESDphysspAveDisc,[minESDphyAve maxESDphyAve])
colorbar_funhan(verticales)
title(hplot,'ESD ave - disc')
grid on
hplot = subplot(2,3,3+3);
imagesc(ESDphysspAveUpperDisc,[minESDphyAveUpper maxESDphyAveUpper])
colorbar_funhan(verticales)
title(hplot,'ESD upper - disc')
grid on
%...................................................................................
%===================================================================================
%...................................................................................
% $$$ jamstecrest_gaussecomodel1D_plotsforLanSmith(logESDphysspAveCont,logESDphysspStdCont,PHYsspcont,ZOOsspcont,DINsspcont,PONsspcont,UXssp,GXssp);
%...................................................................................
%===================================================================================
%...................................................................................
nbins = 4;
%%nbins = 5;
if strcmp(keyModelResol,'0D')
    nbins = ndepths;
end
%...................................................................................
dn = ndepths/nbins;
%...................................................................................
%===================================================================================
%...................................................................................
%%J = [1,[dn:dn:ndepths]];
J = [dn:dn:ndepths];
%...................................................................................
%zzdepths = zdepths + deltaz;
zzdepths = zdepths + deltaz/2; % Adaptation to new vertical levels (Le Gland, 17/09/2019)
zzdepthsJ = zzdepths(J); 
%zzdepthsJ = zzdepthsJ(:); %Column vector.
zzdepthsJ = [0;zzdepthsJ(:)]; % Add zero depth (Le Gland, 17/09/2019)
%...................................................................................
myXtickMarks = [1,[(ndays/12):(ndays/12):ndays]];
%%myXtickMarks = [1,[(ndays/4):(ndays/4):ndays]];
myXtickLabel = num2str(myXtickMarks(:));
%...................................................................................
myYaxisLabel = 'Depth [m]';
%myYtickMarks = [J(:)']; %Add first grid node corresponding to zero meters.
myYtickLabel = [num2str(zzdepthsJ)]; %(for vertical depth resolved 1D model)
%%myYtickMarks = [J(:)']; %Add first grid node corresponding to zero meters.
%%myYtickMarks = [1,J(:)']; %Add first grid node corresponding to zero meters.
myYtickMarks = [0.5,J(:)'+0.5]; %Add first grid node corresponding to zero meters.
%%myYtickLabel = [sprintf('%s 0',0.00);num2str(zzdepthsJ)]; %(for vertical depth resolved 1D model)
%...................................................................................
%===================================================================================
%...................................................................................
PHYmax = max([max(PHYTsspcont(:)),max(PHYTsspdisc(:)),+eps]); 
ZOOmax = max([max(ZOOsspcont(:)),max(ZOOsspdisc(:)),+eps]); 
DINmax = max([max(DINsspcont(:)),max(DINsspdisc(:)),+eps]); 
PONmax = max([max(PONsspcont(:)),max(PONsspdisc(:)),+eps]); 
%...................................................................................
PHYmin = min([min(PHYTsspcont(:)),min(PHYTsspdisc(:)),-eps]); 
ZOOmin = min([min(ZOOsspcont(:)),min(ZOOsspdisc(:)),-eps]); 
DINmin = min([min(DINsspcont(:)),min(DINsspdisc(:)),-eps]); 
PONmin = min([min(PONsspcont(:)),min(PONsspdisc(:)),-eps]); 
%...................................................................................
% $$$ PHYmin = -eps;
% $$$ ZOOmin = -eps;
% $$$ DINmin = -eps;
% $$$ PONmin = -eps;
%...................................................................................
%===================================================================================
%...................................................................................
% $$$ figure(1040)
% $$$ subplot(2,2,1)
% $$$ plot(tode,XAVEode)
% $$$ %%set(gca,'Ylim',[-4 2]) %For x-axis nutrients.
% $$$ xlabel('time [days]')
% $$$ %%title('log(ESD) ave')
% $$$ title('mean')
% $$$ grid on
% $$$ subplot(2,2,2)
% $$$ plot(tode,XSTDode,'b-',tode,XVARode,'r-')
% $$$ %%set(gca,'Ylim',[0 2]) %For x-axis nutrients.
% $$$ xlabel('time [days]')
% $$$ %%title('log(ESD) std -- var')
% $$$ title('std (blue) -- variance (red)')
% $$$ grid on
% $$$ %%return
% $$$ subplot(2,2,3)
% $$$ plot(tode,ESDphysspAveCont)
% $$$ set(gca,'Ylim',[0 10])
% $$$ xlabel('time [days]')
% $$$ title('ESD ave')
% $$$ grid on
% $$$ subplot(2,2,4)
% $$$ plot(tode,ESDphysspAveUpperCont,'b--',tode,ESDphysspAveLowerCont,'b--')
% $$$ set(gca,'Ylim',[0 10])
% $$$ xlabel('time [days]')
% $$$ title('ESD std (Upper / Lower)')
% $$$ grid on
% $$$ print('-dpng','-r300','jamstecrest_gaussecomodel1D_fig002.png')
%...................................................................................
%===================================================================================
% $$$ %...................................................................................
% $$$ logESDaveMax = max(logESDphysspAveCont(:));
% $$$ logESDstdMax = max(logESDphysspStdCont(:));
% $$$ %...................................................................................
% $$$ logESDaveMin = min(logESDphysspAveCont(:));
% $$$ logESDstdMin = min(logESDphysspStdCont(:));
% $$$ %...................................................................................
% $$$ ESDaveMax = max(ESDphysspAveCont(:));
% $$$ ESDstdMax = max(ESDphysspCVcont(:));
% $$$ %...................................................................................
% $$$ ESDaveMin = min(ESDphysspAveCont(:));
% $$$ ESDstdMin = min(ESDphysspCVcont(:));
% $$$ %...................................................................................
%===================================================================================
%...................................................................................
logESDaveMax = max([max(logESDphysspAveCont(:)),max(logESDphysspAveDisc(:))]); 
logESDstdMax = max([max(logESDphysspStdCont(:)),max(logESDphysspStdDisc(:))]); 
%...................................................................................
logESDaveMin = min([min(logESDphysspAveCont(:)),min(logESDphysspAveDisc(:))]); 
logESDstdMin = min([min(logESDphysspStdCont(:)),min(logESDphysspStdDisc(:))]); 
%...................................................................................
ESDaveMax = max([max(ESDphysspAveCont(:)),max(ESDphysspAveDisc(:))]); 
ESDstdMax = max([max(ESDphysspCVcont(:)),max(ESDphysspCVdisc(:))]); 
%...................................................................................
ESDaveMin = min([min(ESDphysspAveCont(:)),min(ESDphysspAveDisc(:))]); 
ESDstdMin = min([min(ESDphysspCVcont(:)),min(ESDphysspCVdisc(:))]); 
%...................................................................................
% Optimal temperature (Le Gland, 02/09/2019)
TOPTaveMax = max([max(TOPTphysspAveCont(:)),max(TOPTphysspAveDisc(:))]); 
TOPTstdMax = max([max(TOPTphysspStdCont(:)),max(TOPTphysspStdDisc(:))]);
%...................................................................................
TOPTaveMin = min([min(TOPTphysspAveCont(:)),min(TOPTphysspAveDisc(:))]); 
TOPTstdMin = min([min(TOPTphysspStdCont(:)),min(TOPTphysspStdDisc(:))]);
%...................................................................................
% Correlation and entropy (Le Gland, 05/11/2019)
CorrelationAbsMax = max([max(abs(physspCorCont(:))),max(abs(physspCorDisc(:)))]);
EntropyMin = min([min(physspEntCont(:)),min(physspEntDisc(:))]);
EntropyMax = max([max(physspEntCont(:)),max(physspEntDisc(:))]);
%===================================================================================
%if strcmp(keyTraitAxis,'ESD')
%    myTitle221 = 'log(ESD) ave';
%    myTitle222 = 'log(ESD) std';
%    myTitle223 = 'ESD ave';
%    %%myTitle224 = 'ESD std';
%    myTitle224 = 'ESD cv';
%elseif strcmp(keyTraitAxis,'SST')
%    myTitle221 = 'SST ave';
%    myTitle222 = 'SST std';
%    myTitle223 = 'exp (SST ave)';
%    myTitle224 = 'exp (SST cv)';
%end
myTitle221 = 'log(ESD) ave';
myTitle222 = 'log(ESD) std';
myTitle223 = 'Topt ave';
myTitle224 = 'Topt std';
myTitle225 = 'Correlation';
myTitle226 = 'Entropy';

%...................................................................................
GALFAssp = galfa*ones(1,ndays);
%NUMUTssp = numut*ones(1,ndays);
%...................................................................................
minSST = min([sstmin,min(isst(:))]);
maxSST = max([sstmax,max(isst(:))]) + eps; %To avoid problems with imagesc limits.
%...................................................................................
minSdin = min(iSdin(:));
maxSdin = max(iSdin(:)) + eps; %To avoid problems with imagesc limits.
%...................................................................................
% $$$ minSdin = min([0,min(iSdin(:))]);
% $$$ maxSdin = max([0,max(iSdin(:))]) + eps; %To avoid problems with imagesc limits.
%...................................................................................
minGalfa = min([1,min(GALFAssp(:))]);
maxGalfa = max([1,max(GALFAssp(:))]) + eps; %To avoid problems with imagesc limits.
%...................................................................................
if strcmp(keyTraitAxis,'ESD')
    myTitle003 = 'DIN supply [mmolN * m-3 * d-1]';
    %%Array2D = iSdin;
    Array2D = SDINssp;
    minArray2D = minSdin;
    maxArray2D = maxSdin;
elseif strcmp(keyTraitAxis,'SST')
    myTitle003 = 'SST variability [Celsius]';
    Array2D = vdepths*isst;
    minArray2D = minSST;
    maxArray2D = maxSST;
end
%...................................................................................
%%UXmax = mup0;
%%GXmax = gzmax;
%...................................................................................
% $$$ UXmax = max(UXssp(:));
% $$$ GXmax = max(GXssp(:));
%...................................................................................
% UXmax = max([max(UXssp(:)),max(uphytotssp(:))]);
% GXmax = max([max(GXssp(:)),max(gphytotssp(:))]);
% Case with 2 traits (Le Gland, 19/07/2019)
UXYmax = max([max(UXYssp(:)),max(uphytotssp(:))]);
GXYmax = max([max(GXYssp(:)),max(gphytotssp(:))]);
%...................................................................................
% $$$ UXmin = min([min(UXssp(:)),min(uphytotssp(:))]);
% $$$ GXmin = min([min(GXssp(:)),min(gphytotssp(:))]);
%...................................................................................
% $$$ UXmin = min(UXssp(:));
% $$$ GXmin = min(GXssp(:));
%...................................................................................
% UXYmin = min(UXYssp(:))-eps;
% GXYmin = min(GXYssp(:))-eps;
% Case with 2 traits (Le Gland, 19/07/2019)
UXYmin = min(UXYssp(:))-eps;
GXYmin = min(GXYssp(:))-eps;
%...................................................................................
% $$$ UXmin = -eps;
% $$$ GXmin = -eps;
%...................................................................................
MUPmax = max([max(MUPsspcont(:)),max(MUPsspdisc(:))]);
MUZmax = max([max(MUZsspcont(:)),max(MUZsspdisc(:))]);
%...................................................................................
MUPmin = min([min(MUPsspcont(:)),min(MUPsspdisc(:))])-eps;
MUZmin = min([min(MUZsspcont(:)),min(MUZsspdisc(:))])-eps;
%...................................................................................
figure(1045)
hplot = subplot(2,2,1);
%%himg = imagesc(NUMUTssp) 
himg = imagesc(GALFAssp,[minGalfa maxGalfa]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'KTW alfa parameter')
grid on
hplot = subplot(2,2,2)
himg = imagesc(log10(PAR2D),log10([0.1 Isat]))
hcbar = colorbar_funhan(verticales);
set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'log10 PAR [W * m-2]')
grid on
%%print('-dpng','-r300','jamstecrest_gaussecomodel1D_fig004.png')
%...................................................................................
%===================================================================================
%%return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT STATISTICS OF BOTH CONTINUOUS AND DISCRETE MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%CONTINOUS:
%...................................................................................
%fignum = 1010;
%...................................................................................
%[hfig1010] = jamstecrest_gaussecomodel1D_imagescuptakerates(MUPsspcont,MUZsspcont,Array2D,NTOTsspcont,fignum,mypackages);
%...................................................................................
%print('-dpng','-r300','jamstecrest_gaussecomodel1D_fig001.png')
%...................................................................................
% Computes the C:Chl ratio in order to ransform biomass into chlorophyll, 
% for comparison with observations, using the algorithm of Lefèvre et al. [2002]
% zangle is aperiodic function with 1 at June solstice and 0 at December solstice
% proxy of the zenithal angle (Le Gland, 15/01/2020)
% Parameters adjusted using Goericke and Welschmeyer. [1998]
zangle = (1/2) * (1 + cos(2*pi*([1:360]-171)/360) );  
% Put the maximum of C:CHL ratio in mid February, to account for a possible delay and for the miwing, not in December (15/07/2020)
%zangle = (1/2) * (1 + cos(2*pi*([1:360]-231)/360) );
Chlratmin = 40;  % mgC/mgChl
Chlratwin = 80; % Surface minimum in winter
Chlratsum = 160; % Surface maximum in summer
Chlrat0   = Chlratwin + (Chlratsum - Chlratwin) * zangle;
Chlratz   = zeros(ndepths,ndays);
for jday = 1:ndays
    Chlratz(:,jday) = Chlratmin + (min(PAR2D(:,jday),25)/25) * (Chlrat0(jday) - Chlratmin);
end
CHLsspcont = PHYTsspcont * (106/16) * 12 ./ Chlratz;
% Model nitracline depth and DIN flux through it (Le Gland, 16/01/2019)
%DINgrad = (DINsspcont(2:end,:) - DINsspcont(1:end-1,:)) / deltaz;
%[DINgradncd,ncdindex] = max(DINgrad,[],1);
%ncd = ncdindex * deltaz;
%DINflux = zeros(1,ndays);
%for jday = 1:ndays
%    DINflux(jday) = DINgradncd(jday) * iKZ(ncdindex(jday)+1,jday);
%end
%DINflux = DINgradncd .* iKZ(ncdindex+1,:);
% We want to avoid "vibration" due to change in nutricline depth and also
% to avoid entrainment effects. Empiracally, we chose the average flux
% between 55m (PP < 55m) and 95m (DIN max > 95m)
DINgrad = (DINsspcont(2:end,:) - DINsspcont(1:end-1,:)) / deltaz;
DINflux = DINgrad .* iKZ(2:end-1,:);
DINsupp = mean(DINflux(6:8,:),1);
%
PPsspcont = PHYTsspcont.*MUPsspcont;
%...................................................................................
fignum = 1020;
%...................................................................................
%[hfig1020] = jamstecrest_gaussecomodel1D_imagescbiomasses(PHYTsspcont,ZOOsspcont,DINsspcont,PONsspcont,fignum,mypackages);
[hfig1020] = jamstecrest_gaussecomodel1D_imagescNPZD(NTOTsspcont,CHLsspcont,PPsspcont,ZOOsspcont,DINsspcont,PONsspcont,fignum,mypackages);
%...................................................................................
% Show observations for comparisons (Le Gland, 11/12/2019)
fignum = 1021;
NTOT_obs = CHL_obs + NO3_obs + PON_obs;
[hfig1021] = jamstecrest_gaussecomodel1D_imagescobs(CHL_obs,PP_obs,NO3_obs,PON_obs,fignum,mypackages);
%%return
%...................................................................................
% Compare model and observations (Le Gland, 07/02/2020)
fignum = 1022;
[hfig1022] = jamstecrest_gaussecomodel1D_imagescmodvsobs(12*(106/16)*PP_obs,CHL_obs,NO3_obs,PON_obs,12*(106/16)*PPsspcont,CHLsspcont,DINsspcont,PONsspcont,fignum,mypackages);
%print('-dpng','-r300','jamstecrest_gaussecomodel1D_fig002.png')
%...................................................................................
fignum = 1030;
%...................................................................................
%[hfig1030] = jamstecrest_gaussecomodel1D_imagescstatistics(logESDphysspAveCont,logESDphysspStdCont,ESDphysspAveCont,ESDphysspCVcont,fignum,mypackages);
% Change figure to plot the two traits, with size only on logscale (Le Gland, 02/09/2019)
[hfig1030] = jamstecrest_gaussecomodel1D_imagescstatistics(logESDphysspAveCont,logESDphysspStdCont,TOPTphysspAveCont,TOPTphysspStdCont,physspCorCont,physspEntCont,fignum,mypackages);
%...................................................................................
%print('-dpng','-r300','jamstecrest_gaussecomodel1D_fig003.png')
%...................................................................................
%===================================================================================
%DISCRETE:
%...................................................................................
%fignum = 2010;
%...................................................................................
%[hfig2010] = jamstecrest_gaussecomodel1D_imagescuptakerates(MUPsspdisc,MUZsspdisc,Array2D,NTOTsspdisc,fignum,mypackages);
%...................................................................................
%print('-dpng','-r300','jamstecrest_discretemodel1D_fig001.png')
%...................................................................................
fignum = 2020;
%...................................................................................
%[hfig2020] = jamstecrest_gaussecomodel1D_imagescbiomasses(PHYTsspdisc,ZOOsspdisc,DINsspdisc,PONsspdisc,fignum,mypackages);
%[hfig2020] = jamstecrest_gaussecomodel1D_imagescNPZD(NTOTsspdisc,PHYTsspdisc,MUPsspdisc,ZOOsspdisc,DINsspdisc,PONsspdisc,fignum,mypackages);
%...................................................................................
%print('-dpng','-r300','jamstecrest_discretemodel1D_fig002.png')
%...................................................................................
fignum = 2030;
%...................................................................................
%[hfig2030] = jamstecrest_gaussecomodel1D_imagescstatistics(logESDphysspAveDisc,logESDphysspStdDisc,ESDphysspAveDisc,ESDphysspCVdisc,fignum,mypackages);
% Change figure to plot the two traits, with size only on logscale (Le Gland, 02/09/2019)
%[hfig2030] = jamstecrest_gaussecomodel1D_imagescstatistics(logESDphysspAveDisc,logESDphysspStdDisc,TOPTphysspAveDisc,TOPTphysspStdDisc,physspCorDisc,physspEntDisc,fignum,mypackages);
%...................................................................................
%print('-dpng','-r300','jamstecrest_discretemodel1D_fig003.png')
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTS FOR OPOS GIJON-IEO AND SWANSEA UNIVERSITY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%SURFACE ONLY VALUES: 
%...................................................................................
% They are the same as xaxis and yaxis ! (Le Gland, 16/09/2019)
xrng = logESDphy;
yrng = ytraitrng; % Le Gland, 06/09/2019
%xrng = ytraitrng; % SST (Le Gland, 23/07/2019) 
%...................................................................................
% $$$ xmin = -4;
% $$$ xmax = +4;
% $$$ xdel = (xmax-xmin)/(nphy-1); 
% $$$ xrng = [xmin:xdel:xmax];
%...................................................................................
jdepth = 1; %Surface only.
%...................................................................................
% phy_cont = ones(ndays,nphy)*nan;
% phy_disc = ones(ndays,nphy)*nan;
% Case with 2 traits (Le Gland, 19/07/2019)
% phy_cont = ones(ndays,nxphy)*nan;
% phy_disc = ones(ndays,nxphy)*nan;
% 2 traits plotted in the same time (Le Gland, 05/09/2019)
phy_cont = ones(ndays,nxphy,nyphy)*nan;
phy_disc = ones(ndays,nxphy,nyphy)*nan;
%...................................................................................
for jday = [1:ndays] 
    phytot = PHYTsspcont(jdepth,jday);
    xave = logESDphysspAveCont(jdepth,jday);
    xsig = logESDphysspStdCont(jdepth,jday);
    % PDF = (1.0 / (xsig * sqrt(2*pi))) * exp( -(xrng - xave).^2 / (2*xsig^2) ); %Okay.
    % phy_cont(jday,:) = phytot .* (PDF .* xdel);
    % 2 traits, SST plotted (Le Gland, 23/07/2019)
    % PDF = (1.0 / (xsig * sqrt(2*pi))) * exp( -(xrng - xave).^2 / (2*xsig^2) ); %Okay.
    % phy_cont(jday,:) = phytot .* (PDF .* xdel);
    % 2D PDF (Le Gland, 05/09/2019)
    yave = TOPTphysspAveCont(jdepth,jday);
    ysig = TOPTphysspStdCont(jdepth,jday);
    xycor = physspCorCont(jdepth,jday);
    %xydet = xsig^2 * ysig^2 - xsig*ysig*xycor;
    xydet = xsig^2 * ysig^2 - (xsig*ysig*xycor)^2; % Le Gland, 19/11/2019
    xymat = [xsig^2, xsig*ysig*xycor; xsig*ysig*xycor, ysig^2];
    xyinvmat = inv(xymat);
    PDF = ones(nxphy,nyphy)*nan;
    for i = 1:nxphy
        for j = 1:nyphy
            % xyvec = [xrng(i) - xave; yrng(j) - yave];
            xyvec = [xaxis(i) - xave; yaxis(j) - yave]; % Use xaxis and yaxis instead of creating new variables (Le Gland, 17/09/2019)
            PDF(i,j) = (1.0 / (sqrt(xydet) * 2 * pi)) * exp( - (1/2) * xyvec' * xyinvmat * xyvec);
        end
    end
    phy_cont(jday,:,:) = phytot .* (PDF .* xdel .* ydel);
end
%...................................................................................
for jday = [1:ndays] 
    phytot = PHYTsspdisc(jdepth,jday);
    xave = logESDphysspAveDisc(jdepth,jday);
    xsig = logESDphysspStdDisc(jdepth,jday);
    PDF = (1.0 / (xsig * sqrt(2*pi))) * exp( -(xrng - xave).^2 / (2*xsig^2) ); %Okay.
    %phy_disc(jday,:) = phytot .* (PDF .* xdel);
    %phy_disc2(jday,:) = PHYsspdisc3D(jdepth,jday,:); % Show real species abundance, not Gaussian equivalent (Le Gland, 16/04/2019)
    %phy_disc2(jday,:) = sum(PHYsspdisc3D(jdepth,jday,:,:),4); % 2 traits (Le Gland, 18/07/2019)
    %phy_disc2(jday,:) = sum(PHYsspdisc3D(jdepth,jday,:,:),3); % SST (Le Gland, 24/07/2019)
    % 2D distribution
    phy_disc(jday,:,:) = PHYsspdisc3D(jdepth,jday,:,:);
end

%...................................................................................
fignum = [70];
%jamstecrest_gaussecomodel1D_samplingplots(phy_cont,phy_disc2,xrng,ndays,fignum);
%jamstecrest_gaussecomodel1D_samplingplots(phy_cont,phy_disc,xrng,yrng,ndays,fignum);
%jamstecrest_gaussecomodel1D_samplingplots(phy_cont,phy_disc,xaxis,yaxis,ndays,fignum);
%jamstecrest_gaussecomodel1D_bellcurve(PHYTsspcont,logESDphysspAveCont,logESDphysspStdCont,TOPTphysspAveCont,TOPTphysspStdCont,physspCorCont,PHYsspdisc3D,xaxis,yaxis,fignum);
jamstecrest_distribution(PHYTsspcont,logESDphysspAveCont,logESDphysspStdCont,TOPTphysspAveCont,TOPTphysspStdCont,physspCorCont,PHYsspdisc3D,xaxis,yaxis,itemp,DINsspdisc,fignum);
%...................................................................................
%print('-dpng ','-r300','jamstecrest_discretemodel1D_plot-phy-lognormal.png')
%print('-depsc','-r300','jamstecrest_discretemodel1D_plot-phy-lognormal.eps')
%...................................................................................
fignum = [80];
jamstecrest_gaussecomodel1D_surftraitplot(DINsspcont,DINsspcont_K,temp(:,1:360),...
PHYTsspcont,PHYTsspcont_K,PHYTsspcont_T,XAVE_sspcont,XAVE_sspcont_K,XSTD_sspcont,XSTD_sspcont_K,...
YAVE_sspcont,YAVE_sspcont_T,YSTD_sspcont,YSTD_sspcont_T,XYCOR_sspcont,myXtickMarks,myXtickLabel,fignum);
%...................................................................................
%===================================================================================
fignum = [24];
jamstecrest_gaussecomodel1D_contvsdiscplot(PHYTsspdisc,PHYTsspcont,logESDphysspAveDisc,logESDphysspAveCont,...
TOPTphysspAveDisc,TOPTphysspAveCont,logESDphysspStdDisc,logESDphysspStdCont,...
TOPTphysspStdDisc,TOPTphysspStdCont,physspCorDisc,physspCorCont,ndepths,fignum);
%...................................................................................
%print('-dpng ','-r300','jamstecrest_gaussecomodel1D_fig024.png')
%...................................................................................
disp('Mean error and root mean square error on total phytoplankton, mean size and standard deviation')
BiasPhy = mean(mean(PHYTsspcont - PHYTsspdisc)) / mean(mean(PHYTsspdisc))
LSQPhy  = sqrt(mean(mean((PHYTsspcont - PHYTsspdisc).^2))) / mean(mean(PHYTsspdisc))
BiasAve = mean(mean(logESDphysspAveCont - logESDphysspAveDisc)) % / mean(mean(logESDphysspAveDisc))
LSQAve  = sqrt(mean(mean((logESDphysspAveCont - logESDphysspAveDisc).^2))) % / mean(mean(logESDphysspAveDisc))
BiasStd = mean(mean(logESDphysspStdCont - logESDphysspStdDisc)) / mean(mean(logESDphysspStdDisc))
LSQStd  = sqrt(mean(mean((logESDphysspStdCont - logESDphysspStdDisc).^2))) / mean(mean(logESDphysspStdDisc))

disp('Total primary production, mean phytoplankton, mean traits, variances and correlation weighted by phytoplankton in the continuous model')
logDIN_minavemax = [min(log(DINsspcont(:))),sum(log(DINsspcont(:)).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(log(DINsspcont(:)))]
T_minavemax = [min(temp(:)),sum(sum(temp(:,1:360).*PHYTsspcont))/sum(PHYTsspcont(:)),max(temp(:))]
PP = deltaz*sum(PHYTsspcont(:).*MUPsspcont(:))*(106/16)*12*0.001 % in gC/m2/yr
Phy_minavemax = [min(PHYTsspcont(:)),mean(PHYTsspcont(:)),max(PHYTsspcont(:))]
logKnave_minavemax = [min(logESDphysspAveCont(:)),sum(logESDphysspAveCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(logESDphysspAveCont(:))]
TOPTave_minavemax = [min(TOPTphysspAveCont(:)),sum(TOPTphysspAveCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(TOPTphysspAveCont(:))]
logKnstd_minavemax = [min(logESDphysspStdCont(:)),sum(logESDphysspStdCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(logESDphysspStdCont(:))]
TOPTstd_minavemax = [min(TOPTphysspStdCont(:)),sum(TOPTphysspStdCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(TOPTphysspStdCont(:))]
Cor_minavemax = [min(physspCorCont(:)),sum(physspCorCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),max(physspCorCont(:))]
%moments = [sum(PHYTsspcont(:).*MUPsspcont(:)),mean(PHYTsspcont(:)),sum(logESDphysspAveCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),...
%           sum(TOPTphysspAveCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),sum(logESDphysspStdCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),...
%           sum(TOPTphysspStdCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:)),sum(physspCorCont(:).*PHYTsspcont(:))/sum(PHYTsspcont(:))]

% Check convergence (Le Gland, 26/02/2020)
%Phy_mov  = Voutcont(Iphy,1:nsteps+1-nsize) -Voutcont(Iphy,nsize:nsteps); % Moving difference with previous year at surface
%xave_mov = Voutcont(Ixave,1:nsteps+1-nsize)./Voutcont(Iphy,1:nsteps+1-nsize) -Voutcont(Ixave,nsize:nsteps);

% Function to check convergence (Le Gland, 17/07/2020)
[conv_P,conv_xave,conv_yave,conv_xxvar,conv_yyvar,conv_xycor] = jamstecrest_convergence_check(PHYTodecont,XAVE_odecont,YAVE_odecont,XXVAR_odecont,YYVAR_odecont,XYCOV_odecont,ndays);

%...................................................................................
return
