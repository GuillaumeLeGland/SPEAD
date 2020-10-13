function [itemp,iparz0,PAR2D,imld,iKZ] = SPEAD_1D_externalForcings(ndepths,zdepths,ndays,kw)

%========================================================================
% Load the external forcings from the input files
%........................................................................
% Temperature is from a monthly climatology (arithmetic mean + interpolation with interp)
% Kz and PAR are from model GOTM (geometric mean)

% Thse forcings are from BATS raw measurements
MLDlevitus360days = load('INPUTS/MLDlevitus360days.txt');
imld = MLDlevitus360days;
PARbats360days = load('INPUTS/PARbats360days.txt');
TEMPbatsinterp2D = load('INPUTS/TEMPbatsinterp2D.txt');

% Interpolate climatological temperature into a ndepths*ndays forcing file
ndepths_temp = size(TEMPbatsinterp2D,1); % Number of depth level in temperature input file (200)
width_temp = floor(ndepths_temp / ndepths); % Ratio of depth levels in the input and forcing files
itemp = zeros(ndepths,ndays);
for i = 1:ndepths
    itemp_month = mean(TEMPbatsinterp2D(width_temp*(i-1)+1:width_temp*i,:),1);
    itemp(i,:) = interp1(1:365,itemp_month,1:364/359:365); % Interpolates on 360 days instead of 365
end

% Load GOTM vertical diffusivity
% ncid  = netcdf.open('INPUTS/phys_BATS.nc','NOWRITE'); % netcdf.open only works in MATLAB
% varid = netcdf.inqVarID(ncid,'nuh');
% Kz1   = netcdf.getVar(ncid,varid);
% netcdf.close(ncid);
% Kz2   = 60*60*24*squeeze(Kz1); % 1x1x250x1095 (m2/s) --> 250x1095 (m2/d)
% Kz3   = flipud(Kz2(51:250,:));
% Kz4   = (Kz3(:,1:365) + Kz3(:,366:730) + Kz3(:,731:1095)) / 3; % Average of three years (200x365) 
% 
% % Interpolate GOTM Kz into a ndepths*ndays forcing file and smooth
% ndepths_Kz = size(Kz4,1); % Number of depth level in temperature input file (200)
% width_Kz   = floor(ndepths_Kz / ndepths); % Ratio of depth levels in the input and forcing files
% iKZ = zeros(ndepths+1,ndays);
% for i = 1:ndepths+1
%     Kz5 = mean(log(Kz4(max(1,width_Kz*(i-3/2)+1):min(width_Kz*(i-1/2),200),:)),1); % Interpolation in log-scale
%     Kz6 = interp1(1:365,Kz5,1:364/359:365); % Interpolates on 360 days instead of 365
%     Kz7 = [Kz6(346:360),Kz6,Kz6(1:15)];
%     Kz8 = smooth(Kz7,30,'sgolay');  % Smoothing
%     iKZ(i,:) = exp(Kz8(16:375)); % Final vertical diffusivity
% end
iKZ = load('INPUTS/KzGOTM.txt');

% Transformation of Einstein/m2/d into W/m2
% Make sure that the iparz0 max is at day 171, and not at day 180 as in PARbats360days
iparz0 = zeros(1,360);
iparz0(1:351)   = 2.5*PARbats360days(10:360);
iparz0(352:360) = 2.5*PARbats360days(1:9);
% iparz0 = 2.5*PARbats360days;

PAR2D = exp(-kw*zdepths(:))*iparz0; % PAR profiles [W*m-2] (depth,time)
%........................................................................
%========================================================================
%************************************************************************
return

