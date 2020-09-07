function [VOLcell,GRMcell,MOLcell]=myConversionCellESDtoMoles(ESD,keyBiomassUnits)

%----------------------------------------------------------------------------
%NOTE:
%Doubling the ESD (equivalent spherical diameter) of a cell should lead to an 
%increase x8 (i.e. x2^3) on the femto (d-15) moles of carbon (MOLcell) per cell.
%----------------------------------------------------------------------------
%Surface Sphere: Surface =   (4)*pi*(Radius)^2
%Volumen Sphere: Volumen = (4/3)*pi*(Radius)^3 = (1/3)*4*pi*(Radius)^3
%Radiuss Sphere: Radiuss = ( Volume/(pi*(4/3)) )^(1/3)
%----------------------------------------------------------------------------
% 1.0 [pgC*um-3] = 1.0 [gC*cm-3]
%----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CARBON DENSITY FACTOR (VOLUMEN TO CARBON MASS CONVERSION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%.......................
fcCarbonToMolles = 1.00/12.00; %[molC*gC-1]
%.......................
%=======================
%.......................
CDFcell = 0.70; %[pgC*um-3] [see BRATBAK 1985, AEM 49(6)]
%.......................
CDFcell = 0.56; %[pgC*um-3] [see BRATBAK 1985, AEM 49(6)]
%.......................
CDFcell = 0.38; %[pgC*um-3] [see LEE and FUHRMAN, AEM, 53(6)]
%.......................
CDFcell = 0.25; %[pgC*um-3] [see Menden-Deuer and Lessard (2000) LIMNOL 45]
%.......................
CDFcell = 0.22; %[pgC*um-3] [see SIERACKI etal. 1995, DSR-I, Fig.3]
%.......................
CDFcell = 0.22; %[pgC*um-3] [see FRY 1988]
%.......................
CDFcell = 0.15; %[pgC*um-3] [see GUNDERSEN etal. 2002, L&0]
%.......................
CDFcell = 0.12; %[pgC*um-3] [see SPITZ etal. 2001, pag.1737]
%.......................
%=======================
%.......................
% $$$ fconv_VolumeToCarbon = 0.15/1.00; %[pgC*um-3] (see Massana96.pdf, Table2)
% $$$ fconv_VolumeToCarbon = 0.50/1.00; %[pgC*um-3] (see Romanova and Sazhin 2010, MarBio)
% $$$ fconv_VolumeToCarbon = 0.56/1.00; %[pgC*um-3] (see BRATBAK 1985, AEM 49(6), pag.Abstract)
% $$$ fcVolumeToCarbon = fconv_VolumeToCarbon*(1.00/10^12); %[gC*um-3]
%.......................
fcVolumeToCarbon = 7.00d-13; %[gC*um-3] [see BRATBAK 1985, AEM 49(6)]
%.......................
%=======================
%.......................
% $$$ ESDcell = 0.25; %[um]
% $$$ ESDcell = 0.30; %[um]
% $$$ ESDcell = 0.50; %[um]
% $$$ ESDcell = 1.00; %[um]
%.......................
% $$$ ESD = [0.25,0.30,0.5,1.0];
%.......................
ESDcell = ESD;
%.......................
%=======================
fconv_mol2fmol = 10^(15)/1; %[mol] x [fmol/mol] = [fmol]
fconv_mol2pmol = 10^(12)/1; %[mol] x [pmol/mol] = [pmol]
%.......................
if strcmp(keyBiomassUnits,'Femto')
    fconvBiomass = fconv_mol2fmol;
elseif strcmp(keyBiomassUnits,'Pico')
    fconvBiomass = fconv_mol2pmol;
end
%.......................
VOLcell = (1/3)*4*pi*(ESDcell/2).^3; %[um3]
GRMcell = VOLcell * fcVolumeToCarbon * fconvBiomass; %[fgC] or [pgC]
MOLcell = GRMcell * fcCarbonToMolles; %[fmolC] or [pmolC]
%.......................
w = [ESDcell;VOLcell;GRMcell;MOLcell]
%.......................
%=======================

return
%*****************************************
ESD = 5.0; %Cell equivalent spherical diameter [um*cell-1]
VOL = (4/3)*pi*(ESD/2)^3; %Cell biovolume [um3*cell-1]
CDF = 0.25;  %Carbon Density Factor [pgC*um-3]
CEL = 1d6; %[cell*m-3]
BIO = CEL*VOL*CDF*(1/1d12)*fcCgramsCmoles*(1000/1); %[mmolC*m-3]
% [cell*m-3] x [um3*cell-1] x [pgC*um-3] x [gC*pgC-1] x [molC*gC] x [mmolC*molC-1] = [m-3] x [mmolC] = [mmolC*m-3] 


fcCgramsCmoles = (1.00/12.00); %[molC*gC-1] 
fcCmolesCmmoles = (1000/1.00); %[mmolC*molC-1]
fcCmolesNmoles = (16.0/106.0); %[molN*molC-1]
fcCpicogramsCgrams = (1/1d12); %[gC*pgC-1]
ESD = 2.0; %Cell equivalent spherical diameter [um*cell-1]
VOL = (4/3)*pi*(ESD/2)^3; %Cell biovolume [um3*cell-1]
CDF = 0.25;  %Carbon Density Factor [pgC*um-3]
%-----------------------------------------------------------------------------------------------------
%[um3*cell-1] x [pgC*um-3] x [gC*pgC-1] x [molC*gC] x [mmolC*molC-1] x [mmolN*molN-1] = [mmolN*cell-1] 
%-----------------------------------------------------------------------------------------------------
BIOcell = (VOL*CDF)*(fcCpicogramsCgrams*fcCgramsCmoles*fcCmolesCmmoles*fcCmolesNmoles); %[mmolN*cell-1]
BIO     = CEL*BIOcell; %[cell*m-3] x [mmolN*cell-1] = [mmolN*m-3]
