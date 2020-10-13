% Script to manage the analysis of model outputs.
% Not the recommanded programming practice, but cleans the main script.
%...................................................................................
%===================================================================================
%...................................................................................
% Phytoplankton growth rate and grazing rate in the 2-trait continuous model
if strcmp(key2T,'yes')
    UXYode = SPEAD_1D_dailyAve(UXYout,deltat);
    GXYode = SPEAD_1D_dailyAve(GXYout,deltat);
end
% 1-trait models
if strcmp(keyKN,'yes')
    UXode = SPEAD_1D_dailyAve(UXout,deltat);
    GXode = SPEAD_1D_dailyAve(GXout,deltat);
end
if strcmp(keyTOPT,'yes')
    UYode = SPEAD_1D_dailyAve(UYout,deltat);
    GYode = SPEAD_1D_dailyAve(GYout,deltat);
end
% Discrete model
if strcmp(keyDisc,'yes')
    uphyodedisc = uphydaydisc; %Phy specific uptake rate of discrete model [d-1]
    gphyodedisc = gphydaydisc; 
    FPHYodedisc = FPHYdaydisc; %Zoo uptake rate of discrete model [mmolN*m-3*d-1]
    GPHYodedisc = GPHYdaydisc; 
end
%...................................................................................
%===================================================================================
%...................................................................................
deltaday = 1;
%...................................................................................
%%Jstep = [1:nsteps]; %For whole simulation (all years at all time steps)
%...................................................................................
Jstep = 1:(ndays/deltat); %For first year only (at all time steps) 
%...................................................................................
%%Jstep = [(ndays/deltat)*(nyear-1)+1:(ndays/deltat)*nyear]; %For last year only (at all time steps) 
%...................................................................................
Jdays = (ndays/deltaday)*(nyear-1)+1:(ndays/deltaday)*nyear; %Last year days (ie. from day 721 until day 1080)
%...................................................................................
%===================================================================================
%...................................................................................
% 2-trait model
if strcmp(key2T,'yes')
    PHYTodecont = Vodecont(Iphy,:);
    ZOOodecont  = Vodecont(Izoo,:);
    DINodecont  = Vodecont(Idin,:);
    PONodecont  = Vodecont(Ipon,:);
    %...............................................................................
    NTOTodecont = PHYTodecont + ZOOodecont + DINodecont + PONodecont;
end
% 1-trait models
if strcmp(keyKN,'yes')
    PHYTodecont_K = Vodecont_K(Iphy,:);
    ZOOodecont_K  = Vodecont_K(Izoo,:);
    DINodecont_K  = Vodecont_K(Idin,:);
    PONodecont_K  = Vodecont_K(Ipon,:);
end
if strcmp(keyTOPT,'yes')
    PHYTodecont_T = Vodecont_T(Iphy,:);
    ZOOodecont_T  = Vodecont_T(Izoo,:);
    DINodecont_T  = Vodecont_T(Idin,:);
    PONodecont_T  = Vodecont_T(Ipon,:);
end
%...................................................................................
%===================================================================================
if strcmp(keyPhysics,'not')
    %...............................................................................
    % Trait distribution moments in the 2-trait continuous model
    if strcmp(key2T,'yes')
        XAVE_odecont  = Vodecont(Ixave,:);
        YAVE_odecont  = Vodecont(Iyave,:);
        XXVAR_odecont = Vodecont(Ixxvar,:);
        YYVAR_odecont = Vodecont(Iyyvar,:);
        XYCOV_odecont = Vodecont(Ixycov,:);
    end
    % 1-trait models
    if strcmp(keyKN,'yes')
        XAVE_odecont_K  = Vodecont_K(Ixave_K,:);
        XXVAR_odecont_K = Vodecont_K(Ixxvar_K,:);
    end
    if strcmp(keyTOPT,'yes')
        YAVE_odecont_T  = Vodecont_T(Iyave_T,:);
        YYVAR_odecont_T = Vodecont_T(Iyyvar_T,:);
    end
    %...............................................................................
elseif strcmp(keyPhysics,'yes')
    %...............................................................................
    % Trait distribution moments in the 2-trait continuous model
    if strcmp(key2T,'yes')
        XAVE_STAR_odecont  = Vodecont(Ixave,:);
        YAVE_STAR_odecont  = Vodecont(Iyave,:);
        XXVAR_STAR_odecont = Vodecont(Ixxvar,:);
        YYVAR_STAR_odecont = Vodecont(Iyyvar,:);
        XYCOV_STAR_odecont = Vodecont(Ixycov,:);
        %...........................................................................
        XAVE_odecont  = (XAVE_STAR_odecont./PHYTodecont);
        YAVE_odecont  = (YAVE_STAR_odecont./PHYTodecont);
        XXVAR_odecont = (XXVAR_STAR_odecont./PHYTodecont) - XAVE_odecont.^2;
        YYVAR_odecont = (YYVAR_STAR_odecont./PHYTodecont) - YAVE_odecont.^2;
        XYCOV_odecont = (XYCOV_STAR_odecont./PHYTodecont) - XAVE_odecont.*YAVE_odecont;
    end
    % 1-trait models
    if strcmp(keyKN,'yes')
        XAVE_STAR_odecont_K  = Vodecont_K(Ixave_K,:);
        XXVAR_STAR_odecont_K = Vodecont_K(Ixxvar_K,:);
        XAVE_odecont_K  = (XAVE_STAR_odecont_K./PHYTodecont_K);
        XXVAR_odecont_K = (XXVAR_STAR_odecont_K./PHYTodecont_K) - XAVE_odecont_K.^2;    
    end
    %...............................................................................
    if strcmp(keyTOPT,'yes')
        YAVE_STAR_odecont_T  = Vodecont_T(Iyave_T,:);
        YYVAR_STAR_odecont_T = Vodecont_T(Iyyvar_T,:);
        YAVE_odecont_T  = (YAVE_STAR_odecont_T./PHYTodecont_T);
        YYVAR_odecont_T = (YYVAR_STAR_odecont_T./PHYTodecont_T) - YAVE_odecont_T.^2;
    end
    %...............................................................................
end
%===================================================================================
%...................................................................................
% Total phytoplankton, zooplankton, DIN and PON in the 2-trait model
if strcmp(key2T,'yes')
    PHYTsspcont = PHYTodecont(:,Jdays);
    ZOOsspcont = ZOOodecont(:,Jdays);
    DINsspcont = DINodecont(:,Jdays);
    PONsspcont = PONodecont(:,Jdays);
    NTOTsspcont = PHYTsspcont + ZOOsspcont + DINsspcont + PONsspcont;
end
% 1-trait models
if strcmp(keyKN,'yes')
    PHYTsspcont_K = PHYTodecont_K(:,Jdays);
    ZOOsspcont_K  = ZOOodecont_K(:,Jdays);
    DINsspcont_K  = DINodecont_K(:,Jdays);
    PONsspcont_K  = PONodecont_K(:,Jdays);
end
if strcmp(keyTOPT,'yes')
    PHYTsspcont_T = PHYTodecont_T(:,Jdays);
    ZOOsspcont_T  = ZOOodecont_T(:,Jdays);
    DINsspcont_T  = DINodecont_T(:,Jdays);
    PONsspcont_T  = PONodecont_T(:,Jdays);
end
%...................................................................................
%===================================================================================
%...................................................................................
% Last year outputs in the 2-trait model
if strcmp(key2T,'yes')
    XAVE_sspcont = XAVE_odecont(:,Jdays);
    YAVE_sspcont  = YAVE_odecont(:,Jdays);
    XXVAR_sspcont = XXVAR_odecont(:,Jdays);
    YYVAR_sspcont = YYVAR_odecont(:,Jdays);
    XYCOV_sspcont = XYCOV_odecont(:,Jdays);
end
% 1-trait models
if strcmp(keyKN,'yes')
    XAVE_sspcont_K  = XAVE_odecont_K(:,Jdays);
    XXVAR_sspcont_K = XXVAR_odecont_K(:,Jdays);
end
if strcmp(keyTOPT,'yes')
    YAVE_sspcont_T  = YAVE_odecont_T(:,Jdays);
    YYVAR_sspcont_T = YYVAR_odecont_T(:,Jdays);
end
%...................................................................................
%===================================================================================
%OBTAIN STANDARD DEVIATION BY SQUARING THE VARIANCE:
%...................................................................................
% 2-trait model
if strcmp(key2T,'yes')
    XSTD_odecont = sqrt(XXVAR_odecont);
    XSTD_sspcont = sqrt(XXVAR_sspcont);
    YSTD_odecont = sqrt(YYVAR_odecont);
    YSTD_sspcont = sqrt(YYVAR_sspcont);
    % Correlation (covariance normalized by standard deviations, always between -1 and 1)
    XYCOR_odecont = XYCOV_odecont ./ (XSTD_odecont .* YSTD_odecont);
    XYCOR_sspcont = XYCOV_sspcont ./ (XSTD_sspcont .* YSTD_sspcont);
end
% 1-trait models
if strcmp(keyKN,'yes')
    XSTD_odecont_K = sqrt(XXVAR_odecont_K);
    XSTD_sspcont_K = sqrt(XXVAR_sspcont_K);
end
if strcmp(keyTOPT,'yes')
    YSTD_odecont_T = sqrt(YYVAR_odecont_T);
    YSTD_sspcont_T = sqrt(YYVAR_sspcont_T);
end
%...................................................................................
%===================================================================================
%CHANGE SOME VARNAMES: 
%...................................................................................
if strcmp(key2T,'yes')
    ESDphysspAveContBis = []; 
    ESDphysspStdContBis = [];
    %...............................................................................
    logESDphysspAveCont = XAVE_sspcont; 
    logESDphysspStdCont = XSTD_sspcont; 
    %...............................................................................
    TOPTphysspAveCont = YAVE_sspcont;
    TOPTphysspStdCont = YSTD_sspcont;
    %...............................................................................
    % Correlations of the continuous model
    physspCorCont = XYCOR_sspcont;
    %...............................................................................
    %===============================================================================
    %...............................................................................
    UXYssp = UXYode(:,Jdays); 
    GXYssp = GXYode(:,Jdays); 
    %...............................................................................
    FPHYTsspcontBis = UXYssp .* PHYTsspcont;
end
%...................................................................................
%===================================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STATE VARIABLES DISCRETE TRAIT MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
if strcmp(keyDisc,'yes')
    [Vodedisc,todedisc] = SPEAD_1D_dailyAve(Voutdisc,deltat);
    %...............................................................................
    %===============================================================================
    %...............................................................................
    PHYodedisc3D = ones(ndepths,ndays*nyear,nxphy,nyphy)*nan;
    %...............................................................................
    PHYodedisc = Vodedisc(Jphy,:); 
    ZOOodedisc = Vodedisc(Jzoo,:);
    DINodedisc = Vodedisc(Jdin,:);
    PONodedisc = Vodedisc(Jpon,:);
    %...............................................................................
    %===============================================================================
    %...............................................................................
    PHYsspdisc3D = ones(ndepths,ndays,nxphy,nyphy)*nan;
    %...............................................................................
    PHYsspdisc = PHYodedisc(:,Jdays);
    ZOOsspdisc = ZOOodedisc(:,Jdays);
    DINsspdisc = DINodedisc(:,Jdays);
    PONsspdisc = PONodedisc(:,Jdays);
    %...............................................................................
    %===============================================================================
    %...............................................................................
    uphysspdisc = uphydaydisc(:,Jdays,:,:);
    gphysspdisc = gphydaydisc(:,Jdays,:,:);
    %...............................................................................
    %===============================================================================
    for iphy = 1:nxphy
        for jphy = 1:nyphy 
            %....................................................................
            Jphyi = Jphy(ndepths*((jphy-1)*nxphy+(iphy-1))+1:ndepths*((jphy-1)*nxphy+iphy));
            %....................................................................
            iPHYssp = PHYsspdisc(Jphyi,:);
            %....................................................................
            PHYsspdisc3D(:,:,iphy,jphy) = iPHYssp; %Phytoplankton [mmolP*m-3] (nphy,depths)
            %....................................................................
            numstrPHYssp = ['PHY',sprintf('%02.0f',iphy),'ssp'];
            %....................................................................
            myassign003 = [numstrPHYssp,' = iPHYssp;'];
            %....................................................................
            eval(myassign003)
            %....................................................................
        end  
    end
    %........................................................................
    PHYTsspdisc = sum(PHYsspdisc3D(:,:,:),3);
    %........................................................................
    FPHYsspdisc3D = uphysspdisc .* PHYsspdisc3D; %[mmolN*m-3*d-1] 
    GPHYsspdisc3D = gphysspdisc .* PHYsspdisc3D; %[mmolN*m-3*d-1] 
    %...............................................................................
    FPHYTsspdisc = sum(FPHYsspdisc3D(:,:,:),3); 
    GPHYTsspdisc = sum(GPHYsspdisc3D(:,:,:),3); 
    %...............................................................................
    uphytotssp = FPHYTsspdisc ./ PHYTsspdisc; %[d-1] 
    gphytotssp = GPHYTsspdisc ./ PHYTsspdisc; %[d-1] 
    %...............................................................................
    NTOTsspdisc = PHYTsspdisc + ZOOsspdisc + DINsspdisc + PONsspdisc;
end
%............................................................................
if strcmp(keyDisc,'yes')
    % Visualize both logESD and Topt
    logESDphyssp3D = ones(ndepths,ndays,nxphy)*nan;
    for iphy = 1:nxphy
        logesdphyi = xaxis(iphy);
        logESDphysspi = logesdphyi*ones(ndepths,ndays);
        logESDphyssp3D(:,:,iphy) = logESDphysspi; %[depth,time,species]
    end 
    TOPTphyssp3D = ones(ndepths,ndays,nyphy)*nan;
    for iphy = 1:nyphy 
        Toptphyi = yaxis(iphy);
        TOPTphysspi = Toptphyi*ones(ndepths,ndays);
        TOPTphyssp3D(:,:,iphy) = TOPTphysspi;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%APPROXIMATE THE CONTINUOUS DISTRIBUTION STATISTICS FROM THE DISCRETE SIMULATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%================================================================================
%................................................................................
if strcmp(keyDisc,'yes')
    [logESDphysspAveDisc,logESDphysspStdDisc,TOPTphysspAveDisc,TOPTphysspStdDisc,physspCorDisc] = SPEAD_1D_covariance(PHYsspdisc3D,logESDphyssp3D,TOPTphyssp3D,nxphy,nyphy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPECIFIC PRODUCTION RATES (PER DAY):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
if strcmp(key2T,'yes')
    MUPsspcont = UXYssp;
    MUZsspcont = GXYssp .* (PHYTsspcont ./ ZOOsspcont);
end
if strcmp(keyDisc,'yes')
    MUPsspdisc = sum(uphysspdisc(:,:,:).*PHYsspdisc3D(:,:,:),3)./PHYTsspdisc;
    MUZsspdisc = sum(gphysspdisc(:,:,:).*PHYsspdisc3D(:,:,:),3)./ZOOsspdisc;
end
%...................................................................................
%===================================================================================
%...................................................................................
nbins = 4;
if strcmp(keyModelResol,'0D')
    nbins = ndepths;
end
%...................................................................................
dn = ndepths/nbins;
%...................................................................................
%===================================================================================
%...................................................................................
J = dn:dn:ndepths;
%...................................................................................
zzdepths = zdepths + deltaz/2;
zzdepthsJ = zzdepths(J); 
zzdepthsJ = [0;zzdepthsJ(:)]; % Add zero depth
%...................................................................................
myXtickMarks = [1,(ndays/12):(ndays/12):ndays];
myXtickLabel = num2str(myXtickMarks(:));
%...................................................................................
myYaxisLabel = 'Depth [m]';
myYtickLabel = num2str(zzdepthsJ); %(for vertical depth resolved 1D model)
myYtickMarks = [0.5,J(:)'+0.5]; %Add first grid node corresponding to zero meters.
%...................................................................................
%===================================================================================
%...................................................................................

% Computes the C:Chl ratio in order to ransform biomass into chlorophyll, 
% for comparison with observations, using the algorithm of Lefèvre et al. [2002]
% zangle is aperiodic function with 1 at June solstice and 0 at December solstice
% proxy of the zenithal angle
% Parameters adjusted using Goericke and Welschmeyer. [1998]
zangle = (1/2) * (1 + cos(2*pi*([1:360]-171)/360) );  
Chlratmin = 40;  % mgC/mgChl
Chlratwin = 80; % Surface minimum in winter
Chlratsum = 160; % Surface maximum in summer
Chlrat0   = Chlratwin + (Chlratsum - Chlratwin) * zangle;
Chlratz   = zeros(ndepths,ndays);
for jday = 1:ndays
    Chlratz(:,jday) = Chlratmin + (min(PAR2D(:,jday),25)/25) * (Chlrat0(jday) - Chlratmin);
end

if strcmp(key2T,'yes')
    CHLsspcont = PHYTsspcont * (106/16) * 12 ./ Chlratz;
    PPsspcont = PHYTsspcont.*MUPsspcont;
end
if strcmp(keyDisc,'yes')
    CHLsspdisc = PHYTsspdisc * (106/16) * 12 ./ Chlratz;
    PPsspdisc = PHYTsspdisc.*MUPsspdisc;
end

if strcmp(keyDisc,'yes') && strcmp(key2T,'yes')
    PHYmax = max([max(PHYTsspcont(:)),max(PHYTsspdisc(:)),+eps]); 
    CHLmax = max([max(CHLsspcont(:)),max(CHLsspdisc(:)),+eps]); 
    ZOOmax = max([max(ZOOsspcont(:)),max(ZOOsspdisc(:)),+eps]); 
    DINmax = max([max(DINsspcont(:)),max(DINsspdisc(:)),+eps]); 
    PONmax = max([max(PONsspcont(:)),max(PONsspdisc(:)),+eps]); 
    %...................................................................................
    PHYmin = min([min(PHYTsspcont(:)),min(PHYTsspdisc(:)),-eps]); 
    CHLmin = min([min(CHLsspcont(:)),min(CHLsspdisc(:)),-eps]); 
    ZOOmin = min([min(ZOOsspcont(:)),min(ZOOsspdisc(:)),-eps]); 
    DINmin = min([min(DINsspcont(:)),min(DINsspdisc(:)),-eps]); 
    PONmin = min([min(PONsspcont(:)),min(PONsspdisc(:)),-eps]);%
    %...................................................................................
    logESDaveMax = max([max(logESDphysspAveCont(:)),max(logESDphysspAveDisc(:)),+eps]);
    logESDaveMin = min([min(logESDphysspAveCont(:)),min(logESDphysspAveDisc(:)),-eps]);
    logESDstdMax = max([max(logESDphysspStdCont(:)),max(logESDphysspStdDisc(:)),+eps]);
    logESDstdMin = min([min(logESDphysspStdCont(:)),min(logESDphysspStdDisc(:)),-eps]);
    TOPTaveMax = max([max(TOPTphysspAveCont(:)),max(TOPTphysspAveDisc(:)),+eps]);
    TOPTaveMin = min([min(TOPTphysspAveCont(:)),min(TOPTphysspAveDisc(:)),-eps]);
    TOPTstdMax = max([max(TOPTphysspStdCont(:)),max(TOPTphysspStdDisc(:)),+eps]);
    TOPTstdMin = min([min(TOPTphysspStdCont(:)),min(TOPTphysspStdDisc(:)),-eps]);
    CorrelationAbsMax = max([max(abs(physspCorCont(:))),max(abs(physspCorDisc(:))),+eps]);
elseif strcmp(keyDisc,'yes')
    PHYmax = max(PHYTsspdisc(:))+eps;
    CHLmax = max(CHLsspdisc(:))+eps;
    ZOOmax = max(ZOOsspdisc(:))+eps; 
    DINmax = max(DINsspdisc(:))+eps; 
    PONmax = max(PONsspdisc(:))+eps; 
    %...................................................................................
    PHYmin = min(PHYTsspdisc(:))-eps;
    CHLmin = min(CHLsspdisc(:))-eps;
    ZOOmin = min(ZOOsspdisc(:))-eps; 
    DINmin = min(DINsspdisc(:))-eps; 
    PONmin = min(PONsspdisc(:))-eps;
    %...................................................................................
    logESDaveMax = max(logESDphysspAveDisc(:))+eps;
    logESDaveMin = min(logESDphysspAveDisc(:))-eps;
    logESDstdMax = max(logESDphysspStdDisc(:))+eps;
    logESDstdMin = min(logESDphysspStdDisc(:))-eps;
    TOPTaveMax = max(TOPTphysspAveDisc(:))+eps;
    TOPTaveMin = min(TOPTphysspAveDisc(:))-eps;
    TOPTstdMax = max(TOPTphysspStdDisc(:))+eps;
    TOPTstdMin = min(TOPTphysspStdDisc(:))-eps;
    CorrelationAbsMax = max(abs(physspCorDisc(:)))+eps;
elseif strcmp(key2T,'yes')
    PHYmax = max(PHYTsspcont(:))+eps; 
    CHLmax = max(CHLsspcont(:))+eps;
    ZOOmax = max(ZOOsspcont(:))+eps; 
    DINmax = max(DINsspcont(:))+eps; 
    PONmax = max(PONsspcont(:))+eps; 
    %...................................................................................
    PHYmin = min(PHYTsspcont(:))-eps;
    CHLmin = min(CHLsspcont(:))-eps; 
    ZOOmin = min(ZOOsspcont(:))-eps; 
    DINmin = min(DINsspcont(:))-eps; 
    PONmin = min(PONsspcont(:))-eps;
    %...................................................................................
    logESDaveMax = max(logESDphysspAveCont(:))+eps;
    logESDaveMin = min(logESDphysspAveCont(:))-eps;
    logESDstdMax = max(logESDphysspStdCont(:))+eps;
    logESDstdMin = min(logESDphysspStdCont(:))-eps;
    TOPTaveMax = max(TOPTphysspAveCont(:))+eps;
    TOPTaveMin = min(TOPTphysspAveCont(:))-eps;
    TOPTstdMax = max(TOPTphysspStdCont(:))+eps;
    TOPTstdMin = min(TOPTphysspStdCont(:))-eps;
    CorrelationAbsMax = max(abs(physspCorCont(:)))+eps;
end
%...................................................................................
GALFAssp = galfa*ones(1,ndays);
%NUMUTssp = numut*ones(1,ndays);
%...................................................................................
minSST = min([ymin,min(itemp(:))]);
maxSST = max([ymax,max(itemp(:))]) + eps; %To avoid problems with imagesc limits.
%...................................................................................
minGalfa = min([1,min(GALFAssp(:))]);
maxGalfa = max([1,max(GALFAssp(:))]) + eps; %To avoid problems with imagesc limits.
%...................................................................................
if strcmp(keyDisc,'yes') && strcmp(key2T,'yes')
    MUPmax = max([max(MUPsspcont(:)),max(MUPsspdisc(:))])+eps;
    MUZmax = max([max(MUZsspcont(:)),max(MUZsspdisc(:))])+eps;
    MUPmin = min([min(MUPsspcont(:)),min(MUPsspdisc(:))])-eps;
    MUZmin = min([min(MUZsspcont(:)),min(MUZsspdisc(:))])-eps;
    %...............................................................................
    UXYmax = max([max(UXYssp(:)),max(uphytotssp(:))])+eps;
    GXYmax = max([max(GXYssp(:)),max(gphytotssp(:))])+eps;
    UXYmax = min([min(UXYssp(:)),min(uphytotssp(:))])+eps;
    GXYmax = min([min(GXYssp(:)),min(gphytotssp(:))])+eps;
elseif strcmp(keyDisc,'yes')
    MUPmax = max(MUPsspdisc(:))+eps;
    MUZmax = max(MUZsspdisc(:))+eps;
    MUPmin = min(MUPsspdisc(:))-eps;
    MUZmin = min(MUZsspdisc(:))-eps;
    %...............................................................................
    UXYmax = max(uphytotssp(:))+eps;
    GXYmax = max(gphytotssp(:))+eps;
    UXYmax = min(uphytotssp(:))+eps;
    GXYmax = min(gphytotssp(:))+eps;
elseif strcmp(key2T,'yes')
    MUPmax = max(MUPsspcont(:))+eps;
    MUZmax = max(MUZsspcont(:))+eps;
    MUPmin = min(MUPsspcont(:))-eps;
    MUZmin = min(MUZsspcont(:))-eps;
    %...............................................................................
    UXYmax = max(UXYssp(:))+eps;
    GXYmax = max(GXYssp(:))+eps;
    UXYmax = min(UXYssp(:))+eps;
    GXYmax = min(GXYssp(:))+eps;
end
%...................................................................................
%===================================================================================
%%return 