function [Xhovmoller,Nhovmoller]=myHovmollerDiagrams(Xarray3D)

%********************************
%Program MYHOVMOLLERDIAGRAMS.m:
%Esta programa obtiene los Hovmoller-Diagrams (pero NO grafica los datos)
%a partir de una matriz 3D de resolucion (1x1)
%********************************
%------------------------------------
%Xarray3D es una matriz 3D de resol (1x1):
%------------------------------------
[msize,nsize,psize]=size(Xarray3D);

% $$$ LATRG = [+90:-10:-90];
% $$$ LATRG = [+90:-05:-90];
LATRG = [+90:-01:-90];
% $$$ LATRG = [+80:-05:-80];

nbins = length(LATRG)-1; %numero de intervalos de latitud (son 18).
PIXELRG = (LATRG-90)*(-1)

[GLOBE,Land] = myLandMaskGlobal;

keyRunMean = 'not';
% $$$ keyRunMean = 'yes';

for jTime = 1:psize
    jTime;
    jG = GLOBE; 
    jX = Xarray3D(:,:,jTime);
    if strcmp(keyRunMean,'not')
	Gk = jG;
	Xk = jX;
	LAT = [1:nbins];
    elseif strcmp(keyRunMean,'yes')
	Gk = [jG(1,:);jG;jG(end,:)]; %Add/repeat first and last row.
	Xk = [jX(1,:);jX;jX(end,:)]; %Add/repeat first and last row.
	LAT = [2:nbins];
    end
    for iLat = LAT(:)'
	if strcmp(keyRunMean,'not')
	    imin = PIXELRG(iLat)+1;
	    imax = PIXELRG(iLat+1);
	elseif strcmp(keyRunMean,'yes')
	    imin = iLat - 1;
	    imax = iLat + 1;
	end
	latminmax = [imin,imax];
	iGband = Gk([imin:imax],:);
	iXband = Xk([imin:imax],:);
	Hsea = find(iGband(:) == 0); %Ocean pixels (i.e. *not* counting land pixels)
	Hpos = find(iXband(:) > 0); %positive numbers.
	narea = length(Hsea);
	nptos = length(Hpos);
	nptosPcnt = nptos/narea;
	if nptosPcnt > 1.0
	    npixels = [narea,nptos,nptosPcnt]
	    disp('Error!!! Pixels with data can *not* more than ocean pixels!!!')
	    pause
	end
	PcntThres = 0.00;
% $$$ 	PcntThres = 0.10;
% $$$ 	PcntThres = 0.20;
	if nptosPcnt > PcntThres
	    iXmean = nanmean(iXband(:));
% $$$ 	    iXmean = nanmedian(iXband(:));
	    iNptos = nptos;
	else 
	    iXmean = nan;
	    iNptos = nan;
	end
	Xhovmoller(iLat,jTime) = iXmean;
	Nhovmoller(iLat,jTime) = iNptos;
    end
end
