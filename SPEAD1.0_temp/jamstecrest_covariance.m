function [XweightAve,XweightStd,YweightAve,YweightStd,XYweightCor] = jamstecrest_covariance(PHY,XTRAIT,YTRAIT)
% Function to plot weighted means, standard deviations and correlations of the traits (Le Gland, 03/09/2019)
global nxphy nyphy
global nsteps ndepths

%===================================================================================
% <http://en.wikipedia.org/wiki/Weighted_arithmetic_mean> 
% <http://en.wikipedia.org/wiki/Mean_square_weighted_deviation> 
%...................................................................................
sumPHY = sum(PHY(:,:,:),3); 
%sumPHYrepmat = repmat(sumPHY,[1,1,nphy]);
sumPHYrepmat = repmat(sumPHY,[1,1,nxphy,nyphy]);
%PHYpcnt = PHY ./ sumPHYrepmat; 
%PHYpcnt = sum(PHY,4) ./ sumPHYrepmat; % 2 traits (Le Gland, 18/07/2019)
%PHYpcnt = squeeze(sum(PHY,3)) ./ sumPHYrepmat; % SST (23/07/2019)
%PHYpcnt = reshape(sum(PHY,3),ndepths,size(PHY,2),nyphy) ./ sumPHYrepmat;
PHYpcnt = PHY ./ sumPHYrepmat;
%...................................................................................
%%Weights = PHY; 
Weights = PHYpcnt; 
%sumWeights = sum(Weights,3);
sumWeights = sum(Weights(:,:,:),3);
%...................................................................................
%COMPUTING AVERAGE OF CELL SIZE WEIGHTED BY BIOMASS CONTRIBUTION: 
% ESDweight = Weights .* ESD; 
%ESDweight = Weights .* sum(ESD,4); % 2 traits (Le Gland, 18/07/2019)
%ESDweight = Weights .* ESD; % SST (23/07/2019)
%ESDweightSum = sum(ESDweight,3); 
%ESDweightAve = ESDweightSum ./ sumWeights; 
%ESDweightAveRepmat = repmat(ESDweightAve,[1,1,nphy]);

XTR = repmat(XTRAIT,[1,1,1,nyphy]);
YTR = permute(repmat(YTRAIT,[1,1,1,nxphy]),[1,2,4,3]);

Xweight = Weights .* XTR;
Yweight = Weights .* YTR;

XweightSum = sum(Xweight(:,:,:),3);
YweightSum = sum(Yweight(:,:,:),3);
XweightAve = XweightSum ./ sumWeights;
YweightAve = YweightSum ./ sumWeights;

% Does not correspond to previous values given by geometricmean ! Why ?)
XweightAveRepmat = repmat(XweightAve,[1,1,nxphy,nyphy]);
Xdispersion = Weights .* (XTR - XweightAveRepmat).^2;
sumXdispersion = sum(Xdispersion(:,:,:),3);
XweightVar = sumXdispersion ./ sumWeights;
XweightStd = sqrt(XweightVar);

YweightAveRepmat = repmat(YweightAve,[1,1,nxphy,nyphy]);
Ydispersion = Weights .* (YTR - YweightAveRepmat).^2;
sumYdispersion = sum(Ydispersion(:,:,:),3);
YweightVar = sumYdispersion ./ sumWeights;
YweightStd = sqrt(YweightVar);

XYdispersion = Weights .* (XTR - XweightAveRepmat) .* (YTR - YweightAveRepmat);
sumXYdispersion = sum(XYdispersion(:,:,:),3);
XYweightCov = sumXYdispersion ./ sumWeights;
XYweightCor = XYweightCov ./ (XweightStd .* YweightStd);

%...................................................................................
%%ESDdispersion = (ESDweight - ESDweightAveRepmat).^2; %WRONG!!!!
%ESDdispersion = Weights .* (ESD - ESDweightAveRepmat).^2; %Okay.
%ESDdispersion = Weights .* (sum(ESD,4) - ESDweightAveRepmat).^2; % 2 traits (Le Gland, 18/07/2019)
%ESDdispersion = Weights .* (ESD - ESDweightAveRepmat).^2; % SST (23/07/2019)
%sumESDdispersion = sum(ESDdispersion,3); 
%...................................................................................
%ESDweightVar = sumESDdispersion ./ sumWeights; 
%ESDweightStd = sqrt(ESDweightVar); 
%...................................................................................
%***********************************************************************************
return




