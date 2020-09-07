function [Hindex,Jindex,Sindex,Eindex,Dindex] = jamstecrest_gaussecomodel1D_ShannonIndex(PHY,mypackages)

%===================================================================================
%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 
%===================================================================================
%...................................................................................
PHY = sum(PHY,4); % 2 traits (Le Gland, 18/07/2019). Also works with 1 trait.
[ndepths,ntimes,nspecies] = size(PHY);
%...................................................................................
Sphymax = nspecies;
Hphymax = log(Sphymax);
%...................................................................................
Sindex = ones(ndepths,ntimes);
Hindex = ones(ndepths,ntimes);
Pindex = ones(ndepths,ntimes,nspecies);
%...................................................................................
%===================================================================================
%FIGURES:
%...................................................................................
for jday = 1:ntimes
    %...............................................................................
    %phy = squeeze(PHY(:,jday,:)); %(depth,species)
    phy = reshape(PHY(:,jday,:),ndepths,nspecies);
    %...............................................................................
    dimsum = 2;
    phytot = sum(phy,dimsum); %Total biomass of all species (depths)
    %...............................................................................
    %phytotrepmat = squeeze(repmat(phytot,[1,1,nspecies]));
    %phytotrepmat = reshape(repmat(phytot,[1,1,nspecies]),ndepths,nspecies)
    phytotrepmat = repmat(phytot,[1,nspecies]);
    %...............................................................................
    probphy = phy./phytotrepmat; %Relative abundance of each species (depths,species)
    %...............................................................................
    LnProbphy = log(probphy); %(depths,species)
    Hj = probphy .* LnProbphy; %(depths,species)
    Hj(find(isnan(Hj))) = 0; %To avoid NANs.
    %...............................................................................
    Ipnct = find(probphy > 0.01); %Where a species represents more than 1 percent.
    OnesZeros = zeros(ndepths,nspecies);
    OnesZeros(Ipnct) = 1.0;
    %...............................................................................
    Sphy = sum(OnesZeros,dimsum); %Species richness (depths) 
    %...............................................................................
    Sindex(:,jday) = Sphy; %Species Richness (> 1pcnt contribution)
    Hindex(:,jday) = -sum(Hj,dimsum); %(depths,time)
    Dindex(:,jday) = 1 ./ sum(probphy.^2,dimsum); %(depths,time)
    Pindex(:,jday,:) = probphy;
    %...............................................................................
end
%...................................................................................
Eindex = exp(Hindex); 
Jindex = Hindex/Hphymax; %[0 - 1]
%...................................................................................
%===================================================================================

%%%%%%%%
%OUTPUT:
%%%%%%%%
%...................................................................................
Hindex = squeeze(Hindex); %Shannon Index.
Jindex = squeeze(Jindex); %Shannon Index normalized.
Sindex = squeeze(Sindex); %Species Richness.
Eindex = squeeze(Eindex); %Exp (Shannon Index).
Dindex = squeeze(Dindex); %Simpson Index.
%...................................................................................
%===================================================================================
%***********************************************************************************
return

figure(501)
subplot(2,2,1)
%%imagesc(Hindex)
imagesc(Dindex)
colorbar_funhan(verticales)
%%title('Shannon index')
title('Simpson index')
grid on
subplot(2,2,2)
imagesc(Jindex)
title('Shannon index normalized [0 - 1]')
colorbar_funhan(verticales)
grid on
subplot(2,2,3)
imagesc(Sindex)
title('Species Richness')
colorbar_funhan(verticales)
grid on
subplot(2,2,4)
imagesc(Eindex)
title('exp (Shannon index)')
colorbar_funhan(verticales)
grid on

maxPHY = max(PHY(:));
figure(1001)
for jphy = 1:nspecies
    subplot(8,8,jphy)
    Phyj = squeeze(PHY(:,:,jphy)); %[depth,time,species]
    %%imagesc(Phyj)
    imagesc(Phyj,[0 maxPHY])
    colorbar_funhan(verticales)
end

figure(2001)
for jphy = 1:nspecies
    subplot(8,8,jphy)
    probPhyj = Pindex(:,:,jphy);
    %%imagesc(probPhyj)
    %%imagesc(probPhyj,[0 1])
    imagesc(probPhyj,[0 0.2])
    colorbar_funhan(verticales)
end
