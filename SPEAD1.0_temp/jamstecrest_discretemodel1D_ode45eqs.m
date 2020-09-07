function [Vdot] = jamstecrest_discretemodel1D_ode45eqs(iTime,V0)
global galfa gbeta
% global gzmax kgz mz betaz betap mpower 
global gzmax kgz mz betaz betap_max betap_xmax betap_xrange mpower
global mp Isat InhFac numutx numuty % Replace numut by 2 different mutation rates (Le Gland, 10/09/2019) 
global alp0 mup0 % knp0 
global amup aknp bmup % aalp
%global Q10 
global Q10a Q10h % Distinct partition coefficients for auto and heterotrophic processes (Le Gland, 31/10/2019)
global mtot0 %discrete model
global temp0 temp % sst Use temperature at all depths (Le Gland, 15/10/2019)
global icounter 
% global todedot 
global keyAssimConstant % keyTraitAxis keyPhysics keySinking % keyAssimConstant added by Le Gland, 03/10/2019
global Sdin Drate % Mdin
global epsPhy omePhy epsZoo omeZoo md 
global deltat % t0 ndays nyear tmax 
global zdepths ndepths deltaz
global KZ % KZI
global parz0
global kw wsink % kp
global Jphy Jzoo Jdin Jpon Jbox %discrete model 
% global nphy nzoo ndin npon Jbox %discrete model 
global nxphy nyphy % nzoo ndin dpon % discrete model with 2 traits (Le Gland, 12/07/2019)
global jday % jjday
global xtrait ytrait 
% global keyFastNumericalSolving
%global keyNutrientSupply 
global keyTraitDiffusion
%global xrng xrngI dx
global dx dy %xrng yrng % 2 traits (Le Gland, 16/07/2019) 
%global KPPXstar KPPXstarI 
%...................................................................................
global uphydaydisc gphydaydisc %OUTPUTS 
%global uphyoutdisc gphyoutdisc %OUTPUTS 
%...................................................................................
global FPHYdaydisc GPHYdaydisc %OUTPUTS 
%global FPHYoutdisc GPHYoutdisc %OUTPUTS 
%...................................................................................
global FPHYToutdisc EPHYToutdisc MPHYToutdisc GPHYToutdisc %OUTPUTS 
global FZOOoutdisc  EZOOoutdisc  MZOOoutdisc 
global FDINoutdisc  FPONoutdisc
%...................................................................................

%%%%%%%%%%%%%%%%%
%STATE VARIABLES:
%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%PHY = ones(ndepths,nphy)*nan;
PHY = ones(ndepths,nxphy,nyphy)*nan;
%PHY = ones(ndepths,nxphy*nyphy)*nan; % Simpler this way ? (Le Gland, 15/07/2019)
%ZOO = ones(ndepths,1)*nan;
%DIN = ones(ndepths,1)*nan;
%PON = ones(ndepths,1)*nan;
%BOX = ones(ndepths,1)*nan;
%...................................................................................
%for iphy = 1:nphy
%  for iphy = 1:nxphy*nyphy   
%     Jphyi = Jphy(ndepths*(iphy-1)+1:ndepths*(iphy));
%     PHYi  = V0(Jphyi);
%     PHY(:,iphy) = PHYi; %Phyplankton [mmolN*m-3] (nphy,depths) 
% end
% Case with 2 traits (Le Gland, 12/07/2019)
% To be checked
% What is Jphy ? Where is it defined ?
for iphy = 1:nxphy
    for jphy = 1:nyphy
        % WRONG ! CAUSE ERROR WHEN nxphy != nyphy
        % Jphyi = Jphy(ndepths*((jphy-1)*nyphy + iphy-1)+1:ndepths*((jphy-1)*nyphy + iphy)); 
        Jphyi = Jphy(ndepths*((jphy-1)*nxphy + iphy-1)+1:ndepths*((jphy-1)*nxphy + iphy));
        PHYi  = V0(Jphyi);
        PHY(:,iphy,jphy) = PHYi;
    end
end
%...................................................................................
ZOO  = V0(Jzoo);
DIN  = V0(Jdin);
PON  = V0(Jpon);
%BOX  = V0(Jbox);
%...................................................................................
%disp(sum(V0(Jphy(:))))
%disp(sum(V0(Jzoo(:))))
%disp(sum(V0(Jdin(:))))
%disp(sum(V0(Jpon(:))))
VNPZD0 = V0([Jphy,Jzoo,Jdin,Jpon]); %Leaving out virtual box.
%...................................................................................
% concs = [PHY,ZOO,DIN,PON];
%concs = [PHY(:,:),ZOO,DIN,PON]; % 2 traits (Le Gland, 15/07/2019)
%disp(size(concs))
%...................................................................................
%===================================================================================
%........................................................................
% PHYT = sum(PHY,2);
PHYT = sum(PHY(:,:),2);

%........................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%
%DAY OF SIMULATION:
%%%%%%%%%%%%%%%%%%%
% Function ot be modified (jday = jjday !) (Le Gland, 15/07/2019)
%[jday,jjday,newday] = jamstecrest_daycounter(jday,iTime);
[jday,newday] = jamstecrest_daycounter(jday,iTime); % jjday is useless (Le Gland, 13/09/2019)

%%%%%%%%%%%
%SHOW TIME:
%%%%%%%%%%%
%===================================================================================
%...................................................................................
% $$$ if mod(iTime,10) == 0 %show every 10 days
% $$$     disp(['runtime: ',num2str(iTime)]) 
% $$$ end
%...................................................................................
% $$$ if mod(iTime,1) == 0 
% $$$     icounter = icounter + 1; 
% $$$ end
%...................................................................................
icounter = floor(iTime/deltat);
%........................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%
%FAST NUMERICAL SOLVING:
%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
% $$$ %........................................................................
% $$$ ZeroValue = 0d0;
% $$$ %........................................................................
% $$$ if strcmp(keyFastNumericalSolving,'yes')
% $$$     %....................................................................
% $$$     Jneg = find( VNPZD0(:) < ZeroValue); 
% $$$     %....................................................................
% $$$     VNPZD0(Jneg) = ZeroValue; %Put zeros where negative concentrations.
% $$$     %....................................................................
% $$$ end
% $$$ %........................................................................
%========================================================================

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
    wconcs = [PHY(:,:),ZOO,DIN,PON] % 2 traits (Le Gland, 16/07/2019) 
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
%disp(mtoti)
%...................................................................................
%%ydistmax = 1d-6; %original. 
ydistmax = 1d-3; %less strict level (for fast-numerical-solving)
ydist = abs(mtoti - mtot0);
%...................................................................................
Inonzero = find(ydist > ydistmax); 
%...................................................................................
%%if strcmp(keyNutrientSupply,'not') 
if length(Inonzero) > 0 
    %masscheck_N = [iTime,jday,jjday,mtot0,mtoti,ydist]
    masscheck_N = [iTime,jday,mtot0,mtoti,ydist]
    disp('Error!!! mass is NOT conserved!')
    pause
end 
%%end 
%...................................................................................
%if mod(jjday,10) == 0 %show every 10 days
if mod(jday,10) == 0 %show every 10 days
    if strcmp(newday,'yes')
    %masscheck_N = [iTime,jday,jjday,mtot0,mtoti,ydist]
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
% FOOD     = ones(ndepths,nphy)*nan; %[depths, phy species]
% FOODalfa = ones(ndepths,nphy)*nan; %[depths, phy species]
% Case with 2 traits (Le Gland, 12/07/2019)
FOOD     = ones(ndepths,nxphy,nyphy)*nan; %[depths, size, Topt]
FOODalfa = ones(ndepths,nxphy,nyphy)*nan;
%...................................................................................
sumFOOD     = ones(ndepths,1)*nan; %[depths]
sumFOODalfa = ones(ndepths,1)*nan; %[depths]
%...................................................................................
FOOD     = PHY;
% FOODalfa = PHY.^(galfa*ones(1,nphy));
%...................................................................................
% sumFOOD     = sum(FOOD,2);
% sumFOODalfa = sum(FOODalfa,2);
% Case with 2 traits (Le Gland, 12/07/2019)
%disp(size(PHY))
%disp(size(galfa))
% DOES NOT WORK, since galfa*ones(1,nxphy,nyphy) is 20*121, and not 20*11*11
%FOODalfa = PHY.^(galfa*ones(1,nxphy,nyphy)); 
%FOODalfa = zeros(ndepths,nxphy,nyphy);
for iphy = 1:nxphy
    for jphy = 1:nyphy
        %FOODalfa(:,iphy,jphy) = PHY(:,iphy,jphy).^galfa;
        % KTW on size only (Le Gland, 27/09/2019)
        FOODalfa(:,iphy,jphy) = PHY(:,iphy,jphy).*((sum(PHY(:,iphy,:),3).^(galfa-1)));
        % FOODalfa(:,iphy,jphy) = PHY(:,iphy,jphy).*((sum(PHY(:,:,jphy),2).^(galfa-1)));
    end
end
sumFOOD     = sum(FOOD(:,:),2);
sumFOODalfa = sum(FOODalfa(:,:),2);
%...................................................................................
% Vmax = (gzmax*ZOO);
% Vmax is temperature-dependent (Le Gland, 04/11/2019)
tempj = temp(:,jday); % (Le Gland, 04/11/2019)
Vmax = (gzmax*ZOO) .* Q10h.^((tempj-temp0)/10);
%...................................................................................
% for iphy = 1:nphy
%     Qswitchi = (FOODalfa(:,iphy) ./ sumFOODalfa); 
%     Qfeeding = (sumFOOD.^gbeta ./ (sumFOOD.^gbeta + kgz^gbeta)); %Okay.
%     Gphy(:,iphy) = Qswitchi .* Qfeeding .* Vmax;
% end 
% Case with 2 traits (Le Gland, 15/07/2019)
Gphy = zeros(ndepths,nxphy,nyphy); % pre-allocate for speed (Le Gland, 26/11/2019)
for iphy = 1:nxphy
    for jphy = 1:nyphy
        Qswitchi = (FOODalfa(:,iphy,jphy) ./ sumFOODalfa);
        Qfeeding = (sumFOOD.^gbeta ./ (sumFOOD.^gbeta + kgz^gbeta)); % Could be out of the loop
        Gphy(:,iphy,jphy) = Qswitchi .* Qfeeding .* Vmax;
    end
end
%...................................................................................
% sumGphy = sum(Gphy,2); %Gross secondary production [mmolN*m-3*d-1]
sumGphy = sum(Gphy(:,:),2); % 2 traits (Le Gland, 15/07/2019)
% Where is Gphy defined ?
%...................................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% $$$ disp('*** pto1 ***')
% $$$ FOODalfa,sumFOODalfa 
% $$$ PHYT = sum(PHY,2); 
% $$$ sumQswitch = sum(Qswitch,2); 
% $$$ sumQswitch,Qfeeding,Vmax
% $$$ sumGphy,PHYT,(sumGphy./PHYT) 
% $$$ pause
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%
%TURBULENT DIFFUSION:
%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%DIFFphy = zeros(ndepths,nphy);
DIFFphy = zeros(ndepths,nxphy,nyphy); % 2 traits (Le Gland, 15/07/2019)
DIFFzoo = zeros(ndepths,1);
DIFFdin = zeros(ndepths,1);
DIFFpon = zeros(ndepths,1);
%...................................................................................
%===================================================================================
%...................................................................................
if ndepths > 1
kz  = KZ (:,jday);
% kzI = KZI(:,jday); % kzI is no longer necessary (Le Gland, 10/09/2019)
%...................................................................................
% for iphy = 1:nphy
%     %...............................................................................
%     [DIFFphyi] = jamstecrest_TurbulentDiffusion(PHY(:,iphy),deltat,deltaz,kz,kzI,ndepths,'Reflectante');
%     %...............................................................................
%     DIFFphy(:,iphy) = DIFFphyi;
%     %...............................................................................
% end
% 2 traits (Le Gland, 15/07/2019)
for iphy = 1:nxphy
    for jphy = 1:nyphy
        %...............................................................................
        % [DIFFphyi] = jamstecrest_TurbulentDiffusion(PHY(:,iphy,jphy),deltat,deltaz,kz,kzI,ndepths,'Reflectante');
        % kzI is no longer required (Le Gland, 10/09/2019)
        [DIFFphyi] = jamstecrest_TurbulentDiffusion(PHY(:,iphy,jphy),deltat,deltaz,kz,ndepths,'Reflectante');
        %...............................................................................
        DIFFphy(:,iphy,jphy) = DIFFphyi;
        %...............................................................................
    end
end
%...................................................................................
% [DIFFzoo] = jamstecrest_TurbulentDiffusion(ZOO,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
% [DIFFdin] = jamstecrest_TurbulentDiffusion(DIN,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
% [DIFFpon] = jamstecrest_TurbulentDiffusion(PON,deltat,deltaz,kz,kzI,ndepths,'Reflectante');
% kzI is no longer required (Le Gland, 10/09/2019)
[DIFFzoo] = jamstecrest_TurbulentDiffusion(ZOO,deltat,deltaz,kz,ndepths,'Reflectante');
[DIFFdin] = jamstecrest_TurbulentDiffusion(DIN,deltat,deltaz,kz,ndepths,'Reflectante');
[DIFFpon] = jamstecrest_TurbulentDiffusion(PON,deltat,deltaz,kz,ndepths,'Reflectante');
end
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%
%VERTICAL SINKING:
%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
% ADVphy = zeros(ndepths,1);
% ADVzoo = zeros(ndepths,1);
% ADVdin = zeros(ndepths,1);
% ADVpon = zeros(ndepths,1);
%........................................................................
%========================================================================
%........................................................................
[ADVpon] = jamstecrest_SinkingAdvection(PON,deltaz,wsink,ndepths);
% Transform PON to DIN at bottom to avoid PON accumulation (Le Gland, 11/12/2019)
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
%........................................................................
% $$$ SSTj = SST(jday,:); 
% $$$ %........................................................................
% $$$ Q10sst = Q10 .^ ((SSTj - sst0) / 10); %perfil de sst limitation para dia-i.
%........................................................................
% Q10sst = 1.0;
%Q10sst = Q10 ^((sst(jday) - sst0)/10); % Le Gland (05/06/2019), to be improved
%Q10sst = Q10.^((temp(:,jday) - sst0)/10); % Le Gland (15/10/2019), to be improved
%........................................................................
%========================================================================
%........................................................................
% $$$ alp = alp0*ones(ndepths,1); 
% $$$ mup = mup0*ones(ndepths,1);
%........................................................................
% $$$ alp = alp0*Qpar; 
% $$$ mup = mup0*Qpar;
%........................................................................
% Q10 is now based on Topt, not on in-situ temperature (Le Gland, 23/10/2019)
alp = alp0.*Qpar; % .*Q10sst; 
mup = mup0.*Qpar; % .*Q10sst;
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHYTOPLANKTON NUTRIENT UPTAKE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for iphy = 1:nphy 

% if strcmp(keyTraitAxis,'ESD')
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %PHYTOPLANKTON GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %===============================================================================
%     %...............................................................................
%     knp = (mup./alp); 
%     %...............................................................................
%     logESDphyi = xtrait(iphy); %log[um]
%     % logknpphyi = xtrait(iphy); % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     %...............................................................................
%     %===============================================================================
%     %EXPONENTIAL UPTAKE AFFINITY: 
%     %...............................................................................
%     knphy = knp .* exp(aknp*logESDphyi); %Phy half-sat uptake as a function of cell size [mmolN*m-3]
%     % knphy = exp(logknpphyi); % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     %...............................................................................
%     %===============================================================================
%     %EXPONENTIAL MAXIMUM GROWTH RATE: 
%     %...............................................................................
% % $$$     muphy = mup .* exp(amup*logESDphyi); %Phy maximum grazing rate as a function of cell size [d-1]
%     %...............................................................................
%     %===============================================================================
%     %UNIMODAL MAXIMUM GROWTH RATE: 
%     %...............................................................................
%     muphy = mup .* exp(amup*logESDphyi + bmup*logESDphyi.^2);
%     % muphy = mup .* exp(amup*(logknpphyi-log(knp))/aknp); % Case where trait is half-saturation (Le Gland, 26/06/2019)
%     %...............................................................................
%     %===============================================================================
%     %...............................................................................
%     Qdim = (DIN  ./ (DIN + knphy)); %[n.d.]
%     %...............................................................................
%     Uphy(:,iphy) = muphy .* Qdim; %Uptake rate [d-1]  
%     %...............................................................................
%     %===============================================================================
% 
% elseif strcmp(keyTraitAxis,'SST')
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %USING NEW TRAIT Y-AXIS FOR THE OPTIMAL TEMPERATURE FOR GROWTH:
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %===============================================================================
%     %...............................................................................
%     knp = knp0*ones(ndepths,1);
%     %...............................................................................
%     %%Toptphyi = ytrait(iphy); %[Celsius]
%     Toptphyi = xtrait(iphy); % Le Gland, 05/06/2019
%     %...............................................................................
%     gammay = 4.0; %Tolerance range [Celsius]
%     %...............................................................................
%     sstj = sst(jday); % sstj has to be defined first ! (Le Gland, 25/04/2019) 
%     Qsst = exp(-(Toptphyi - sstj).^2 / (2*gammay^2)); %SST niches [0 - 1]
%     %...............................................................................
%     Qdim = (DIN  ./ (DIN + knp)); %[n.d.]
%     %...............................................................................
%     Uphy(:,iphy) = mup .* Qdim .* Qsst; %Uptake rate [d-1] 
%     %...............................................................................
%     %===============================================================================
% end %endif

% %...................................................................................
% Fphy(:,iphy) = Uphy(:,iphy) .* PHY(:,iphy); %Phy production [mmolN*m-3*d-1] 
% % Test with grazing equal to growth (Le Gland, 09/05/2019) 
% % Gphy(:,iphy) = Fphy(:,iphy);
% %...................................................................................
% %%Mphy(:,iphy) =   mp  * PHY(:,iphy); %Phy mortality [mmolN*m-3*d-1] %Linear.
% %...................................................................................
% Mphy(:,iphy) =   mp  * (PHY(:,iphy) .* PHYT); %Quadratic. %Okay.
% % Test with grazing equal to growth (Le Gland, 09/05/2019)
% % Mphy(:,iphy) = 0;
% %...................................................................................
% %%Mphy(:,iphy) =   mp  * (PHYT.^mpower); %Quadratic. %WRONG!!!
% %...................................................................................
% end %endfor iphy = 1:nphy

% Case with 2 traits (Le Gland, 15/07/2019)
betap = zeros(nxphy,nyphy); % Assimilation efficiency for every species (Le Gland, 11/11/2019)
% Pre-allocate for efficiency (Le Gland, 26/11/2019)
Uphy = zeros(ndepths,nxphy,nyphy);
Fphy = zeros(ndepths,nxphy,nyphy);
Mphy = zeros(ndepths,nxphy,nyphy);
Ephy = zeros(ndepths,nxphy,nyphy);
for iphy = 1:nxphy
    for jphy = 1:nyphy

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PHYTOPLANKTON GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %===============================================================================
        %...............................................................................
        knp = (mup./alp); 
        %...............................................................................
        logESDphyi = xtrait(iphy); %log[um]
        Toptphyi   = ytrait(jphy);
        %===============================================================================
        %EXPONENTIAL UPTAKE AFFINITY: 
        %...............................................................................
        knphy = knp .* exp(aknp*logESDphyi); %Phy half-sat uptake as a function of cell size [mmolN*m-3]
        %...............................................................................
        %===============================================================================
        %UNIMODAL MAXIMUM GROWTH RATE: 
        %...............................................................................
        muphy = mup .* exp(amup*logESDphyi + bmup*logESDphyi.^2);
        % muphy = mup .* exp(amup*(logknpphyi-log(knp))/aknp); % Case where trait is half-saturation (Le Gland, 26/06/2019)
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
        %gammay = 4.0; %Tolerance range [Celsius]
        %...............................................................................
        %sstj = sst(jday); % sstj has to be defined first ! (Le Gland, 25/04/2019)
        %sstj = temp(:,jday); % Use temperature at all depths (Le Gland, 15/10/2019)
        % Qsst = exp(-(Toptphyi - sstj).^2 / (2*gammay^2)); %SST niches [0 - 1]
        % New skewed response to temperature (Le Gland, 23/10/2019)
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
        % Test with grazing equal to growth (Le Gland, 09/05/2019) 
        % Gphy(:,iphy) = Fphy(:,iphy);
        %...................................................................................
        % Mphy(:,iphy,jphy) =   mp  * PHY(:,iphy,jphy); %Phy mortality [mmolN*m-3*d-1] %Linear.
        Mphy(:,iphy,jphy) =   mp  * PHY(:,iphy,jphy) .* Q10h.^((tempj-temp0)/10); % Temperature dependence on mortality (Le Gland, 11/11/2019)
        % Mphy(:,iphy,jphy) =   mp  * PHY(:,iphy,jphy) .* PHYT(:) .* Q10h.^((tempj-temp0)/10); % Quadratic
        %...................................................................................
        % Mphy(:,iphy,jphy) =   mp  * (PHY(:,iphy,jphy) .* PHYT); %Quadratic. %Okay.
        
        if strcmp(keyAssimConstant,'not')
            betap(iphy,jphy) = betap_max .* exp( -1/2 .* (logESDphyi - betap_xmax).^2 ./ (betap_xrange.^2) );
        else
            betap(iphy,jphy) = betap_max;
        end
        Ephy(:,iphy,jphy) = (1 - betap(iphy,jphy)) .* Fphy(:,iphy,jphy);
        
    end
end
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ZOOPLANKTON GRAZING AND PRODUCTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
% Fzoo = sum(Gphy,2); %Gross secondary production [mmolN*m-3*d-1]
Fzoo = sum(Gphy(:,:),2); % 2 traits (Le Gland, 15/07/2019)
Ezoo = (1-betaz) * Fzoo; %Zooplankton exudation [mmolN*m-3*d-1]
%Ephy = (1-betap) * Fphy; %Phy exudation [mmolN*m-3*d-1] 
% Test with grazing equal to growth (Le Gland, 09/05/2019)
% Ephy = 0 * Fphy;
%Mzoo = mz*(ZOO.^mpower); %Zoo mortality [mmolN*m-3*d-1]
Mzoo = mz*(ZOO.^mpower) .* Q10h.^((tempj-temp0)/10); % Temperature dependence on mortality (Le Gland, 11/11/2019)

%...................................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% $$$ disp('*** pto2 ***')
% $$$ Fzoo./sum(PHY,2) 
% $$$ sum(PHY,2) 
% $$$ Fzoo 
% $$$ pause 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%===================================================================================
%...................................................................................
% Mpon = md*PON;
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
% sumFphy = sum(Fphy,2); 
% sumEphy = sum(Ephy,2); 
% sumMphy = sum(Mphy,2);
% Case with 2 traits (Le Gland, 15/07/2019)
sumFphy = sum(Fphy(:,:),2); 
sumEphy = sum(Ephy(:,:),2); 
sumMphy = sum(Mphy(:,:),2);

%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%
%TRAIT DIFFUSION:
%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%Ndiff = zeros(ndepths,nphy);
Ndiff = zeros(ndepths,nxphy,nyphy); % 2 traits (Le Gland, 17/07/2019)
%...................................................................................
if strcmp(keyTraitDiffusion,'yes')
    for jdepth = 1:ndepths
        %...................................................................................
        % N = PHY(jdepth,:);
        % rx = Uphy(jdepth,:);
        % Case with 2 traits (Le Gland, 15/07/2019)
        N    = squeeze(PHY(jdepth,:,:));
        % rxy  = squeeze(Uphy(jdepth,:,:));
        % Trait diffusion has to take exudation into acccount, since
        % exudation limits reproduction (Le Gland, 04/11/2019)
        % rxy = squeeze(Uphy(jdepth,:,:) - Ephy(jdepth,:,:)./PHY(jdepth,:,:));
        rxy = squeeze(Uphy(jdepth,:,:)) .* betap(:,:);
        % nuxy = numut(jdepth);
        nux = numutx(jdepth);
        nuy = numuty(jdepth);
        %...................................................................................
        % $$$     %%BC_type = 'Absorbente';
        % $$$     BC_type = 'Reflectante';
        % $$$     NI = interp1(xrng,N,xrngI);
        % $$$     kppx  = KPPXstar(jdepth,:) .*N; %[m2*d-1]
        % $$$     kppxI = KPPXstarI(jdepth,:).*NI; %[m2*d-1]
        % $$$     [jNdiff] = jamstecrest_TurbulentDiffusion(N,deltat,dx,kppx,kppxI,nphy,BC_type)';
        %...................................................................................
        %[jNdiff] = jamstecrest_discretemodel1D_traitDispersion(N,rx,nux,nphy,dx);
        %[jNdiff] = jamstecrest_discretemodel1D_traitDispersion(N,rxy,nuxy,nxphy,nyphy,dx,dy);
        % One mutation rate for each trait (Le Gland, 10/09/2019)
        [jNdiff] = jamstecrest_discretemodel1D_traitDispersion(N,rxy,nux,nuy,nxphy,nyphy,dx,dy);
        %...................................................................................
        %Ndiff(jdepth,:) = jNdiff; %size(ndepths,nphy): Trait-diffusion across phyto size-class at each vertical depth.
        Ndiff(jdepth,:,:) = jNdiff; % 2 traits (Le Gland, 17/07/2019)
        %...................................................................................
    end
end
%...................................................................................
%===================================================================================


%===================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ORDINARY DIFFERENTIAL EQUATIONS (ODEs):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
% $$$ %...................................................................................
% $$$ for iphy = 1:nphy 
% $$$     dPHYdt(:,iphy) = Fphy(:,iphy) - Gphy(:,iphy) - Mphy(:,iphy) + Ndiff(:,iphy); 
% $$$ end
% $$$ %...................................................................................
% $$$ dZOOdt =      Fzoo - Ezoo - Mzoo; 
% $$$ %...................................................................................
% $$$ dDINdt = - sumFphy + sumMphy*   omePhy  + Mzoo*   omeZoo  + Ezoo*   epsZoo  + Mpon;
% $$$ %...................................................................................
% $$$ dPONdt =           + sumMphy*(1-omePhy) + Mzoo*(1-omeZoo) + Ezoo*(1-epsZoo) - Mpon;
% $$$ %...................................................................................
%===================================================================================
%...................................................................................
% for iphy = 1:nphy 
%     dPHYdt(:,iphy) = Fphy(:,iphy) - Ephy(:,iphy) - Gphy(:,iphy) - Mphy(:,iphy) + Ndiff(:,iphy); 
% end
% Case with 2 traits (Le Gland, 16/07/2019)
for iphy = 1:nxphy
    for jphy = 1:nyphy 
        dPHYdt(:,iphy,jphy) = Fphy(:,iphy,jphy) - Ephy(:,iphy,jphy) - Gphy(:,iphy,jphy) - Mphy(:,iphy,jphy) + Ndiff(:,iphy,jphy); 
    end
end
%...................................................................................
dZOOdt =   Fzoo - Ezoo - Mzoo; 
%...................................................................................
dDINdt = - sumFphy + epsPhy*sumEphy +     omePhy*sumMphy +     epsZoo*Ezoo +     omeZoo*Mzoo + Mpon; 
%...................................................................................
dPONdt =         (1-epsPhy)*sumEphy + (1-omePhy)*sumMphy + (1-epsZoo)*Ezoo + (1-omeZoo)*Mzoo - Mpon; 
%...................................................................................
%===================================================================================
%...................................................................................
FPHYToutdisc(:,icounter) = sumFphy;
EPHYToutdisc(:,icounter) = sumEphy;
GPHYToutdisc(:,icounter) = sumGphy;
MPHYToutdisc(:,icounter) = sumMphy;
%...................................................................................
FZOOoutdisc(:,icounter) = Fzoo;
EZOOoutdisc(:,icounter) = Ezoo;
MZOOoutdisc(:,icounter) = Mzoo;
%...................................................................................
FDINoutdisc(:,icounter) = epsPhy*sumEphy + omePhy*sumMphy + epsZoo*Ezoo + omeZoo*Mzoo + Mpon; 
FPONoutdisc(:,icounter) = (1-epsPhy)*sumEphy + (1-omePhy)*sumMphy + (1-epsZoo)*Ezoo + (1-omeZoo)*Mzoo;
%...................................................................................
%===================================================================================
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%DEBUGGING:
% $$$ %...................................................................................
% $$$ % $$$ PHY,Uphy,Gphy,Fphy,Mphy
% $$$ % $$$ pause
% $$$ % $$$ %...................................................................................
% $$$ % $$$ dWdt = [dPHYdt,dZOOdt,dDINdt,dPONdt]
% $$$ % $$$ pause
% $$$ %...................................................................................
% $$$ PHYT = sum(PHY,2); 
% $$$ dPHYTdt = sum(dPHYdt,2); 
% $$$ sumEphy = zeros(ndepths,1); 
% $$$ %...................................................................................
% $$$ masscheck_N = [iTime,jday,mtot0,mtoti,ydist]
% $$$ wGdx = [   Fzoo./PHYT] %[d-1] 
% $$$ wUdx = [sumFphy./PHYT] %[d-1] 
% $$$ wBIOs = [PHYT,ZOO,DIN,PON] 
% $$$ wODEs = [dPHYTdt,dZOOdt,dDINdt,dPONdt] 
% $$$ %...................................................................................
% $$$ dispPHYdot = [dPHYTdt,sumFphy,sumEphy,Fzoo,sumMphy]
% $$$ dispZOOdot = [dZOOdt,Fzoo,Ezoo,Mzoo]
% $$$ dispDINdot = [dDINdt,sumFphy,epsPhy*sumEphy,omePhy*sumMphy,epsZoo*Ezoo,omeZoo*Mzoo,Mpon] 
% $$$ dispPONdot = [dPONdt,(1-epsPhy)*sumEphy,(1-omePhy)*sumMphy,(1-epsZoo)*Ezoo,(1-omeZoo)*Mzoo,Mpon] 
% $$$ %...................................................................................
% $$$ disp('-------------')
% $$$ pause
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK FLUX CONSERVATION:
%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%sumODEs = sum(dPHYdt,2) + dZOOdt + dDINdt + dPONdt;
sumODEs = sum(dPHYdt(:,:),2) + dZOOdt + dDINdt + dPONdt; % 2 traits (Le Gland, 16/07/2019)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADD EXTERNAL NUTRIENT SOURCES AND ALL STATE VARIABLES SINKS BY DILUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
% $$$ %...................................................................................
% $$$ Dphy = Drate(jday)*PHY; %Dilution of phytoplankton [mmolN*m-3*d-1]
% $$$ Dzoo = Drate(jday)*ZOO;
% $$$ Ddin = Drate(jday)*DIN;
% $$$ Dpon = Drate(jday)*PON;
% $$$ %...................................................................................
% $$$ Sno3 = Sdin(jday)*ones(ndepths,1);
% $$$ %...................................................................................
%===================================================================================
%...................................................................................
Dphy = Drate(icounter)*PHY; %Dilution of phytoplankton [mmolN*m-3*d-1]
%disp(Dphy(1,7,11))
Dzoo = Drate(icounter)*ZOO;
Ddin = Drate(icounter)*DIN;
Dpon = Drate(icounter)*PON;
%...................................................................................
Sno3 = Sdin(:,icounter); %Supply of nutrients into the system [mmolN*m-3*d-1]
%...................................................................................
%%Dzoo = zeros(ndepths,1); %Remove washout dilution of zooplanton.
%...................................................................................
% sumDphy = sum(Dphy,2); 
sumDphy = sum(Dphy(:,:),2); % 2 traits (Le Gland, 16/07/2019)
%...................................................................................
%===================================================================================
%...................................................................................
dPHYdt = dPHYdt - Dphy;
dZOOdt = dZOOdt - Dzoo;
dDINdt = dDINdt - Ddin + Sno3; %Only DIN has an input flux source.
dPONdt = dPONdt - Dpon;
%...................................................................................
dBOXdt = -Sno3 + sumDphy + Dzoo + Ddin + Dpon; %Virtual box to check mass conservation.
%...................................................................................
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
%...................................................................................
BOXdot = dBOXdt; 
%...................................................................................
%===================================================================================

%%%%%%%%%%
%STOCKAGE:
%%%%%%%%%%
%===================================================================================
%...................................................................................
% $$$ todedotday(1,jday) = jday;
%...................................................................................
% Where is uphydaydisc defined ? Risk to create a bug with 2 traits (Le Gland, 16/07/2019)
% uphydaydisc is only defined as a "global" in jamstecrest_gaussecomodel1D
% Its dimensions are not specified.
%uphydaydisc(:,jday,:) = Fphy./PHY; %size(depth,time,species) [d-1]
%gphydaydisc(:,jday,:) = Gphy./PHY;
%%...................................................................................
%FPHYdaydisc(:,jday,:) = Fphy;
%GPHYdaydisc(:,jday,:) = Gphy;
% Case with 2 traits (Le Gland, 17/07/2019)
uphydaydisc(:,jday,:,:) = Fphy./PHY; %size(depth,time,species) [d-1]
gphydaydisc(:,jday,:,:) = Gphy./PHY;
%...................................................................................
FPHYdaydisc(:,jday,:,:) = Fphy;
GPHYdaydisc(:,jday,:,:) = Gphy;

%...................................................................................
%===================================================================================
%TOO SLOW TO RUN...
% $$$ %...................................................................................
% $$$ % $$$ todedotoutdisc(1,icounter) = icounter*deltat;
% $$$ %...................................................................................
% $$$ uphyoutdisc(:,icounter,:) = Fphy./PHY; %size(depth,time,species) [d-1]
% $$$ gphyoutdisc(:,icounter,:) = Gphy./PHY;
% $$$ %...................................................................................
% $$$ FPHYoutdisc(:,icounter,:) = Fphy;
% $$$ GPHYoutdisc(:,icounter,:) = Gphy;
% $$$ %...................................................................................
%===================================================================================

%%%%%%%%
%OUTPUT:
%%%%%%%%
%...................................................................................
Vdot = [PHYdot(:);ZOOdot;DINdot;PONdot;BOXdot];
%...................................................................................
%***********************************************************************************
return

