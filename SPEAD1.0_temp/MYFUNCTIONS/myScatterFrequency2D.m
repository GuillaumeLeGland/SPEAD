function [NPTOS] = myScatterFrequency2D(xdata,ydata,phy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MY METHOD PICO-Y-PALA (FASTER):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=================================================
%.................................................
nbins = 128;
mbins = 128;
%.................................................
%=================================================
%.................................................
XI = xdata(:)';
YI = ydata(:)';
%.................................................
%=================================================
%.................................................
msize = numel(XI);
%.................................................
Xmin = min(xdata);
Xmax = max(xdata);
%..................................................
Ymin = min(ydata);
Ymax = max(ydata);
%..................................................
%=================================================
%..................................................
XYmin = [Xmin,Ymin]
XYmax = [Xmax,Ymax]
%..................................................
XL = linspace(Xmin,Xmax,mbins);
YL = linspace(Ymin,Ymax,nbins);
%..................................................
FXL = [1:mbins];
FYL = [1:nbins];
%..................................................
FXI = interp1(XL,FXL,XI,'nearest');
FYI = interp1(YL,FYL,YI,'nearest');
%..................................................
%=================================================
%..................................................
XCOOR = ones(nbins,mbins)*nan;
YCOOR = ones(nbins,mbins)*nan;
%..................................................
NPTOS = ones(nbins,mbins)*nan;
%..................................................
phy_ave = ones(msize,1)*nan;
fxy_frequency = ones(msize,1)*nan;
%..................................................
Xbins = XL;
Ybins = YL;
%..................................................
for jbin = 1:mbins
    jbindisplay = jbin 
    xcoor_j = Xbins(jbin);
    for ibin = 1:nbins
	ycoor_i = Ybins(ibin);
	I = find(FXI == jbin & FYI == ibin);
	nptos = length(I);
	if nptos > 0
	    fxy_frequency(I) = nptos;
	    phy_ave(I) = median(phy(I));
	    NPTOS(ibin,jbin) = nptos;
	end
	XCOOR(ibin,jbin) = xcoor_j;
	YCOOR(ibin,jbin) = ycoor_i;
    end
end
%..................................................
%=================================================
%*************************************************
return

%..................................................
NPTOS = flipud(NPTOS); %To use with "imagesc" (comment out if using "pcolor")
XCOOR = flipud(XCOOR);
YCOOR = flipud(YCOOR);
%..................................................
NPTOS_star = 100*(NPTOS/max(NPTOS(:))); %Normalized percentage [0 - 1]*100
%..................................................

cmap = jet;
figure(30)
subplot(2,2,2)
scatter(xdata,ydata,10,log10(fxy_frequency),'filled')
%%fscatter(xdata,ydata,log10(fxy_frequency),cmap);
colorbar('vertic')


S = phy_ave/max(phy_ave(:))*100;
figure(40)
subplot(2,2,2)
scatter(xdata,ydata,S,phy_ave,'filled')
%%fscatter(xdata,ydata,phy_ave,cmap);
colorbar('vertic')





