function [ESDweightAve,ESDweightStd] = jamstecrest_geometricmean(PHY,ESD)
%global nphy 
global nxphy nyphy % 2 traits (Le Gland, 18/07/2019)
global nsteps ndepths

%===================================================================================
% <http://en.wikipedia.org/wiki/Weighted_arithmetic_mean> 
% <http://en.wikipedia.org/wiki/Mean_square_weighted_deviation> 
%...................................................................................
nphy = nyphy; % 2 traits (Le Gland, 18/07/2019)
%sumPHY = sum(PHY,3);
sumPHY = sum(PHY(:,:,:),3); 
sumPHYrepmat = repmat(sumPHY,[1,1,nphy]); 
%PHYpcnt = PHY ./ sumPHYrepmat; 
%PHYpcnt = sum(PHY,4) ./ sumPHYrepmat; % 2 traits (Le Gland, 18/07/2019)
%PHYpcnt = squeeze(sum(PHY,3)) ./ sumPHYrepmat; % SST (23/07/2019)
PHYpcnt = reshape(sum(PHY,3),ndepths,size(PHY,2),nyphy) ./ sumPHYrepmat;
%...................................................................................
%%Weights = PHY; 
Weights = PHYpcnt; 
sumWeights = sum(Weights,3); 
%...................................................................................
%COMPUTING AVERAGE OF CELL SIZE WEIGHTED BY BIOMASS CONTRIBUTION: 
% ESDweight = Weights .* ESD; 
%ESDweight = Weights .* sum(ESD,4); % 2 traits (Le Gland, 18/07/2019)
ESDweight = Weights .* ESD; % SST (23/07/2019)
ESDweightSum = sum(ESDweight,3); 
ESDweightAve = ESDweightSum ./ sumWeights; 
ESDweightAveRepmat = repmat(ESDweightAve,[1,1,nphy]); 
%...................................................................................
%%ESDdispersion = (ESDweight - ESDweightAveRepmat).^2; %WRONG!!!!
%ESDdispersion = Weights .* (ESD - ESDweightAveRepmat).^2; %Okay.
%ESDdispersion = Weights .* (sum(ESD,4) - ESDweightAveRepmat).^2; % 2 traits (Le Gland, 18/07/2019)
ESDdispersion = Weights .* (ESD - ESDweightAveRepmat).^2; % SST (23/07/2019)
sumESDdispersion = sum(ESDdispersion,3); 
%...................................................................................
ESDweightVar = sumESDdispersion ./ sumWeights; 
ESDweightStd = sqrt(ESDweightVar); 
%...................................................................................
%===================================================================================
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%...................................................................................
% $$$ %%ESDweightAveBis = mean(ESDweight,3); %WRONG!!!!
% $$$ %%ESDweightAveBis = wmean(ESD,Weights,3); %Okay.
% $$$ %%ESDweightStdBis = std(ESDweight,[],3); %WRONG!!!!
% $$$ %%ESDweightStdBis = wstd(Weights,ESD);
% $$$ %%[ESDweightAveBis,ESDweightStdBis] = wmeanstd(ESD,Weights,3); 
% $$$ [ESDweightAveBis,ESDweightStdBis] = wstdmean(ESD,Weights,3,'Biased'); %Okay.
%...................................................................................
% $$$ figure(10)
% $$$ subplot(2,2,1)
% $$$ imagesc(ESDweightAve)
% $$$ colorbar 
% $$$ grid on
% $$$ subplot(2,2,2)
% $$$ imagesc(ESDweightStd)
% $$$ colorbar 
% $$$ grid on
% $$$ subplot(2,2,3)
% $$$ imagesc(ESDweightAveBis)
% $$$ colorbar 
% $$$ grid on
% $$$ subplot(2,2,4)
% $$$ imagesc(ESDweightStdBis)
% $$$ colorbar 
% $$$ grid on
%...................................................................................
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%===================================================================================
%***********************************************************************************
return




