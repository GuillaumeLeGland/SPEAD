function [Hindex,Jindex,Sindex,Eindex,Dindex] = myShannonIndex(PHY) %PHY(180,360,12,64)
%***************************************************************************************
%Use: [Hindex,Jindex,Sindex,Eindex] = myShannonIndex(PHY) %PHY(lat,lon,time,species)
%***************************************************************************************

%%[mlat,nlon,ptime,qspecies] = size(PHY);

mlat = 20;
nlon =  1;
ptime = 360;
qspecies = 64;

Sphymax = qspecies;
Hphymax = log(Sphymax);

Sindex = ones(mlat,nlon,ptime);
Hindex = ones(mlat,nlon,ptime);
Pindex = ones(mlat,nlon,ptime,qspecies);
for k=1:ptime
      monthk = k
      if nlon > 1
	  PHYk = squeeze(PHY(:,:,k,:)); %(lat,long,species)
	  dimsum = 3;
      elseif nlon == 1
	  PHYk = squeeze(PHY(:,k,:)); %(depth,species)
	  dimsum = 2;
      end
      sumPHYk = sum(PHYk,dimsum); %Total biomass of all species (180,360)
      sumPHYkrepmat = squeeze(repmat(sumPHYk,[1,1,qspecies]));

      probPHYk = PHYk./sumPHYkrepmat; %Relative abundance of each species (180,360,64)
      Ipnct = find(probPHYk > 0.01); %Where a species represents more than 1 percent.

      W = squeeze(zeros(mlat,nlon,qspecies));
      W(Ipnct) = 1.0;

      Sphyk = sum(W,dimsum); %Species richness (180,360) 
      LnProbPHYk = log(probPHYk); %(180,360,64)
      Hk = probPHYk .* LnProbPHYk; %(180,360,64)
      Hk(find(isnan(Hk))) = 0; %To avoid NANs.

      Sindex(:,:,k) = Sphyk;
      Hindex(:,:,k) = -sum(Hk,dimsum); %(180,360)
      Dindex(:,:,k) = 1 ./ sum(probPHYk.^2,dimsum); %(180,360)
      Pindex(:,:,k,:) = probPHYk;
end
Eindex = exp(Hindex); 
Jindex = Hindex/Hphymax; %[0 - 1]

%%%%%%%%
%OUTPUT:
%%%%%%%%
Hindex = squeeze(Hindex); %Shannon Index.
Jindex = squeeze(Jindex); %Shannon Index normalized.
Sindex = squeeze(Sindex); %Species Richness.
Eindex = squeeze(Eindex); %Exp (Shannon Index).
Dindex = squeeze(Dindex); %Simpson Index.
%*********************************
return

figure(500)
subplot(2,2,1)
%%imagesc(Hindex)
imagesc(Dindex)
mycolorbar
%%title('Shannon index')
title('Simpson index')
grid on
subplot(2,2,2)
imagesc(Jindex)
title('Shannon index normalized [0 - 1]')
mycolorbar
grid on
subplot(2,2,3)
imagesc(Sindex)
title('Species Richness')
mycolorbar
grid on
subplot(2,2,4)
imagesc(Eindex)
title('exp (Shannon index)')
mycolorbar
grid on

maxPHY = max(PHY(:));
figure(1000)
for jphy = 1:qspecies
    subplot(8,8,jphy)
    Phyj = squeeze(PHY(:,:,jphy)); %[depth,time,species]
    %%imagesc(Phyj)
    imagesc(Phyj,[0 maxPHY])
    mycolorbar
end

figure(2000)
for jphy = 1:qspecies
    subplot(8,8,jphy)
    probPhyj = squeeze(Pindex(:,:,:,jphy));
    %%imagesc(probPhyj)
    %%imagesc(probPhyj,[0 1])
    imagesc(probPhyj,[0 0.2])
    mycolorbar
end
