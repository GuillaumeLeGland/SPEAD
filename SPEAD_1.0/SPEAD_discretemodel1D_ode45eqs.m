function [Vdot] = SPEAD_discretemodel1D_ode45eqs(iTime,V0)
global galfa gbeta
global gzmax kgz mz betaz mpower
global mp Isat InhFac numutx numuty 
global alp0 mup0 
global amup aknp
global Q10a Q10h % Distinct partition coefficients for auto and heterotrophic processes
global mtot0 %discrete model total mass
global temp0 temp
global icounter 
global omePhy epsZoo omeZoo md 
global deltat
global zdepths ndepths deltaz
global KZ
global parz0
global kw wsink
global Jphy Jzoo Jdin Jpon 
global nxphy nyphy
global jday
global xaxis yaxis
global keyTraitDiffusion
global xdel ydel 
%...................................................................................
global uphydaydisc gphydaydisc %OUTPUTS 
%...................................................................................
global FPHYdaydisc GPHYdaydisc %OUTPUTS 
%...................................................................................
global FPHYToutdisc MPHYToutdisc GPHYToutdisc %OUTPUTS 
global FZOOoutdisc  EZOOoutdisc  MZOOoutdisc 
global FDINoutdisc  FPONoutdisc
%...................................................................................

%%%%%%%%%%%%%%%%%
%STATE VARIABLES:
%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
PHY = ones(ndepths,nxphy,nyphy)*nan;
for iphy = 1:nxphy
    for jphy = 1:nyphy
        Jphyi = Jphy(ndepths*((jphy-1)*nxphy + iphy-1)+1:ndepths*((jphy-1)*nxphy + iphy));
        PHYi  = V0(Jphyi);
        PHY(:,iphy,jphy) = PHYi;
    end
end
%...................................................................................
ZOO  = V0(Jzoo);
DIN  = V0(Jdin);
PON  = V0(Jpon);
%...................................................................................
VNPZD0 = V0([Jphy,Jzoo,Jdin,Jpon]);
%...................................................................................
%===================================================================================
%........................................................................
% PHYT = sum(PHY(:,:),2);
%........................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%
%DAY OF SIMULATION:
%%%%%%%%%%%%%%%%%%%
[jday,newday] = SPEAD_1D_daycounter(jday,iTime);

%%%%%%%%%%%
%SHOW TIME:
%%%%%%%%%%%
%===================================================================================
%...................................................................................
icounter = floor(iTime/deltat);
%........................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK FOR NEGATIVE CONCENTRATIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
Ineg = find(VNPZD0 < 0); 
%........................................................................
if length(Ineg > 0)
    pause
    iTime
    disp(['P , Z , N , D']);
    %wconcs = [PHY,ZOO,DIN,PON]
    wconcs = [PHY(:,:),ZOO,DIN,PON] 
    wconcsNeg = VNPZD0(Ineg); 
    disp('Error!!! there are NEGATIVE concentrations!')
    pause
end
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK MASS CONSERVATION:
%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
mtoti = sum(V0); %Checking mass conservation.
%...................................................................................
%%ydistmax = 1d-6;
ydistmax = 1d-3; %less strict level (for fast-numerical-solving)
ydist = abs(mtoti - mtot0);
%...................................................................................
Inonzero = find(ydist > ydistmax); 
%...................................................................................
if length(Inonzero) > 0 
    masscheck_N = [iTime,jday,mtot0,mtoti,ydist]
    disp('Error!!! mass is NOT conserved!')
    pause
end 
%...................................................................................
if mod(jday,10) == 0 %show every 10 days
    if strcmp(newday,'yes')
    masscheck_N = [iTime,jday,mtot0,mtoti,ydist]
    display('-------------------------------------------------')
    end
end
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GRAZING FUNCTIONAL RESPONSE KTW:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
% FOOD     = ones(ndepths,nxphy,nyphy)*nan; %[depths, size, Topt]
FOODalfa = ones(ndepths,nxphy,nyphy)*nan;
%...................................................................................
sumFOOD     = ones(ndepths,1)*nan; %[depths]
% sumFOODalfa = ones(ndepths,1)*nan; %[depths]
%...................................................................................
% FOOD     = PHY;
for iphy = 1:nxphy
    for jphy = 1:nyphy
        FOODalfa(:,iphy,jphy) = PHY(:,iphy,jphy).^galfa;
        % KTW on size only
        %FOODalfa(:,iphy,jphy) = PHY(:,iphy,jphy).*((sum(PHY(:,iphy,:),3).^(galfa-1)));
        % FOODalfa(:,iphy,jphy) = PHY(:,iphy,jphy).*((sum(PHY(:,:,jphy),2).^(galfa-1)));
    end
end
sumFOOD     = sum(PHY(:,:),2);
sumFOODalfa = sum(FOODalfa(:,:),2);
%...................................................................................
% Vmax is temperature-dependent
tempj = temp(:,jday);
Vmax = (gzmax*ZOO) .* Q10h.^((tempj-temp0)/10);
%...................................................................................
Gphy = zeros(ndepths,nxphy,nyphy); % pre-allocate for speed
for iphy = 1:nxphy
    for jphy = 1:nyphy
        Qswitchi = (FOODalfa(:,iphy,jphy) ./ sumFOODalfa);
        Qfeeding = (sumFOOD.^gbeta ./ (sumFOOD.^gbeta + kgz^gbeta)); % Could be out of the loop
        Gphy(:,iphy,jphy) = Qswitchi .* Qfeeding .* Vmax;
    end
end
%................................................................................... 
sumGphy = sum(Gphy(:,:),2); %Gross secondary production [mmolN*m-3*d-1]
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%
%TURBULENT DIFFUSION:
%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
DIFFphy = zeros(ndepths,nxphy,nyphy);
DIFFzoo = zeros(ndepths,1);
DIFFdin = zeros(ndepths,1);
DIFFpon = zeros(ndepths,1);
%...................................................................................
%===================================================================================
%...................................................................................
if ndepths > 1
kz  = KZ (:,jday);
for iphy = 1:nxphy
    for jphy = 1:nyphy
        %...............................................................................
        [DIFFphyi] = SPEAD_1D_TurbulentDiffusion(PHY(:,iphy,jphy),deltat,deltaz,kz,ndepths,'Implicit',100);
        %...............................................................................
        DIFFphy(:,iphy,jphy) = DIFFphyi;
        %...............................................................................
    end
end
%...................................................................................
[DIFFzoo] = SPEAD_1D_TurbulentDiffusion(ZOO,deltat,deltaz,kz,ndepths,'Implicit',100);
[DIFFdin] = SPEAD_1D_TurbulentDiffusion(DIN,deltat,deltaz,kz,ndepths,'Implicit',100);
[DIFFpon] = SPEAD_1D_TurbulentDiffusion(PON,deltat,deltaz,kz,ndepths,'Implicit',100);
end
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%
%VERTICAL SINKING:
%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
[ADVpon] = SPEAD_1D_SinkingAdvection(PON,deltaz,wsink,ndepths);
% Transform PON to DIN at bottom to avoid PON accumulation
ADVdin = zeros(ndepths,1);
ADVdin(end) = (wsink/deltaz)*PON(end);
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
% Use normalized Follows (2007) formula, with photoinhibition of large cells
Kpar   = log(InhFac+1) / Isat;
Kinhib = Kpar / InhFac; 
Fmax   = (Kpar+Kinhib)/Kpar * exp( -(Kinhib/Kpar) * log(Kinhib/(Kpar+Kinhib)) );  
% Qpar is normalized to have a maximum of 1 at Isat
Qpar = Fmax * (1 - exp(-Kpar*jPAR)) .* exp(-Kinhib*jPAR);
%........................................................................
%========================================================================
%........................................................................
%% alp = alp0*ones(ndepths,1); 
%% mup = mup0*ones(ndepths,1);
%........................................................................
alp = alp0.*Qpar;
mup = mup0.*Qpar;
%........................................................................
%========================================================================

% Pre-allocate for efficiency
Uphy = zeros(ndepths,nxphy,nyphy);
Fphy = zeros(ndepths,nxphy,nyphy);
Mphy = zeros(ndepths,nxphy,nyphy);
for iphy = 1:nxphy
    for jphy = 1:nyphy

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PHYTOPLANKTON GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %===============================================================================
        %...............................................................................
        knp = (mup./alp); 
        %...............................................................................
        logESDphyi = xaxis(iphy); %log[um]
        Toptphyi   = yaxis(jphy);
        %===============================================================================
        %EXPONENTIAL UPTAKE AFFINITY: 
        %...............................................................................
        knphy = knp .* exp(aknp*logESDphyi); %Phy half-sat uptake as a function of cell size [mmolN*m-3]
        %...............................................................................
        %===============================================================================
        %UNIMODAL MAXIMUM GROWTH RATE: 
        %...............................................................................
        muphy = mup .* exp(amup*logESDphyi);
        %...............................................................................
        %===============================================================================
        %...............................................................................
        Qdim = (DIN  ./ (DIN + knphy)); %[n.d.] 
        %...............................................................................
        %===============================================================================
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %USING NEW TRAIT Y-AXIS FOR THE OPTIMAL TEMPERATURE FOR GROWTH:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %===============================================================================
        %...............................................................................
        if tempj < Toptphyi + 5
            Qsst = exp(0.2*(tempj-Toptphyi)) .* Q10a.^((Toptphyi-temp0)/10) .* (Toptphyi + 5 - tempj)/5;
        else 
            Qsst = 0;
        end
        % Qsst = exp(-(Toptphyi - tempj).^2 / (2*gammay^2)) .* Q10a.^((Toptphyi-temp0)/10); % Symmetrical Qsst (Le Gland, 15/11/2019)
        %...............................................................................
        Uphy(:,iphy,jphy) = muphy .* Qdim .* Qsst; % Uptake rate [d-1]
        %...............................................................................
        %===============================================================================

        %...................................................................................
        Fphy(:,iphy,jphy) = Uphy(:,iphy,jphy) .* PHY(:,iphy,jphy); %Phy production [mmolN*m-3*d-1] 
        %...................................................................................
        % Mphy(:,iphy,jphy) =   mp  *  PHY(:,iphy,jphy); %Phy mortality [mmolN*m-3*d-1] %Linear -- Okay
        % Mphy(:,iphy,jphy) =   mp  * (PHY(:,iphy,jphy) .* PHYT); %Quadratic -- Okay
        % Mphy(:,iphy,jphy) =   mp  *  PHY(:,iphy,jphy) .* PHYT(:) .* Q10h.^((tempj-temp0)/10); % Quadratic
          Mphy(:,iphy,jphy) =   mp  *  PHY(:,iphy,jphy) .* Q10h.^((tempj-temp0)/10); % Temperature dependence on mortality
        %...................................................................................
    end
end
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ZOOPLANKTON GRAZING AND PRODUCTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
Fzoo = sum(Gphy(:,:),2);
Ezoo = (1-betaz) * Fzoo;
Mzoo = mz*(ZOO.^mpower) .* Q10h.^((tempj-temp0)/10); % Temperature dependence on mortality
%...................................................................................
%===================================================================================
%...................................................................................
% Temperature-dependent, based on heterotrophic Q10 (Le Gland, 04/11/2019)
mdt = md*Q10h.^((tempj-temp0)/10);
Mpon = mdt.*PON;
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TOTAL PHY PRIMARY PRODUCTION, EXUDATION AND MORTALITY: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
sumFphy = sum(Fphy(:,:),2);  
sumMphy = sum(Mphy(:,:),2);
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%
%TRAIT DIFFUSION:
%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
Ndiff = zeros(ndepths,nxphy,nyphy); % 2 traits (Le Gland, 17/07/2019)
%...................................................................................
if strcmp(keyTraitDiffusion,'yes')
    for jdepth = 1:ndepths
        %...................................................................................
        N    = squeeze(PHY(jdepth,:,:));
        rxy = squeeze(Uphy(jdepth,:,:));
        nux = numutx(jdepth);
        nuy = numuty(jdepth);
        %...................................................................................
        [jNdiff] = SPEAD_discretemodel1D_traitDispersion(N,rxy,nux,nuy,nxphy,nyphy,xdel,ydel);
        %...................................................................................
        Ndiff(jdepth,:,:) = jNdiff;
        %Ndiff(jdepth,:,:) = 0.0;
        %...................................................................................
    end
end
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ORDINARY DIFFERENTIAL EQUATIONS (ODEs):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
for iphy = 1:nxphy
    for jphy = 1:nyphy 
        dPHYdt(:,iphy,jphy) = Fphy(:,iphy,jphy) - Gphy(:,iphy,jphy) - Mphy(:,iphy,jphy) + Ndiff(:,iphy,jphy); 
    end
end
%...................................................................................
dZOOdt =   Fzoo - Ezoo - Mzoo; 
%...................................................................................
dDINdt = - sumFphy +     omePhy*sumMphy +     epsZoo*Ezoo +     omeZoo*Mzoo + Mpon; 
%...................................................................................
dPONdt =             (1-omePhy)*sumMphy + (1-epsZoo)*Ezoo + (1-omeZoo)*Mzoo - Mpon; 
%...................................................................................
dBOXdt = sum(sum(dPHYdt,3),2) + dZOOdt + dDINdt + dPONdt; %Virtual box to check mass conservation.
%...................................................................................
%===================================================================================
%...................................................................................
FPHYToutdisc(:,icounter) = sumFphy;
GPHYToutdisc(:,icounter) = sumGphy;
MPHYToutdisc(:,icounter) = sumMphy;
%...................................................................................
FZOOoutdisc(:,icounter) = Fzoo;
EZOOoutdisc(:,icounter) = Ezoo;
MZOOoutdisc(:,icounter) = Mzoo;
%...................................................................................
FDINoutdisc(:,icounter) = omePhy*sumMphy + epsZoo*Ezoo + omeZoo*Mzoo + Mpon; 
FPONoutdisc(:,icounter) = (1-omePhy)*sumMphy + (1-epsZoo)*Ezoo + (1-omeZoo)*Mzoo;
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK FLUX CONSERVATION:
%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
sumODEs = sum(dPHYdt(:,:),2) + dZOOdt + dDINdt + dPONdt;
%...................................................................................
xdistmax = 1d-6;
xdist = abs(sumODEs - 0d0);
%...................................................................................
Jnonzero = find(xdist > xdistmax);
%...................................................................................
if length(Jnonzero) > 1
    display(jcnp)
    sumODEs %the sum of all ODEs should be equal to zero (closed system)!
    disp('Error!!!! Zero flux is NOT conserved!!!')
    pause
end
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%
%ADD PHYSICAL PROCESSES:
%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
PHYdot = dPHYdt + DIFFphy;
ZOOdot = dZOOdt + DIFFzoo;
DINdot = dDINdt + DIFFdin + ADVdin; % Bottom remineralization (ADVdin) added (Le Gland, 11/12/2019)
PONdot = dPONdt + DIFFpon + ADVpon; %Only PON has vertical sinking.
BOXdot = dBOXdt; %Virtual box (only for mass conversation checking)
%...................................................................................
%===================================================================================

%%%%%%%%%%
%STOCKAGE:
%%%%%%%%%%
%===================================================================================
%...................................................................................
uphydaydisc(:,jday,:,:) = Fphy./PHY; %size(depth,time,species) [d-1]
gphydaydisc(:,jday,:,:) = Gphy./PHY;
%...................................................................................
FPHYdaydisc(:,jday,:,:) = Fphy;
GPHYdaydisc(:,jday,:,:) = Gphy;
%...................................................................................
%===================================================================================

%%%%%%%%
%OUTPUT:
%%%%%%%%
%...................................................................................
Vdot = [PHYdot(:);ZOOdot;DINdot;PONdot;BOXdot];
%...................................................................................
%***********************************************************************************
return
