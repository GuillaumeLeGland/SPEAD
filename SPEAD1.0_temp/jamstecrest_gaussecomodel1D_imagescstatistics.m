%function [hfig] = jamstecrest_gaussecomodel1D_imagescstatistics(XAVEssp,XSTDssp,ESDave,ESDcv,fignum,mypackages)
%function [hfig] = jamstecrest_gaussecomodel1D_imagescstatistics(XAVEssp,XSTDssp,YAVEssp,YSTDssp,XYCORssp,XYENTssp,fignum,mypackages)
function [hfig] = jamstecrest_gaussecomodel1D_imagescstatistics(XAVEssp,XSTDssp,YAVEssp,YSTDssp,XYCORssp,PHYssp,fignum,mypackages)
%...................................................................................
global myTitle221 myTitle222 myTitle223 myTitle224  myTitle225 myTitle226
global myXtickMarks myXtickLabel myXaxisLabel
global myYtickMarks myYtickLabel myYaxisLabel
global zdepths ndepths deltaz
global t0 deltat ndays nyear tmax 
%...................................................................................
global logESDaveMax logESDstdMax 
global logESDaveMin logESDstdMin 
global    ESDaveMax    ESDstdMax 
global    ESDaveMin    ESDstdMin 
global   TOPTaveMax   TOPTstdMax % Le Gland, 02/09/2019
global   TOPTaveMin   TOPTstdMin % Le Gland, 02/09/2019
global         CorrelationAbsMax % Le Gland, 05/11/2019 
global   EntropyMin   EntropyMax % Le Gland, 05/11/2019
%...................................................................................

%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 

%FIGURES:
hfig = figure(fignum);
%hplot = subplot(2,2,1);
hplot = subplot(2,3,2);
%himg = imagesc(XAVEssp,[logESDaveMin,logESDaveMax]);
%himg = imagesc(XAVEssp,[-0.8,+0.6]);
himg = imagesc(XAVEssp,[min(XAVEssp(:)),max(XAVEssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
% Do not cut first and last row in half (Le Gland, 17/09/2019) 
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
%title(hplot,myTitle221)
title(hplot,'log(Kn) ave')
grid on

%hplot = subplot(2,2,2);
hplot = subplot(2,3,5);
%himg = imagesc(XSTDssp,[logESDstdMin,logESDstdMax]);
himg = imagesc(XSTDssp,[min(XSTDssp(:)),max(XSTDssp(:))]);
%himg = imagesc(XSTDssp,[0.3,0.6]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
% Do not cut first and last row in half (Le Gland, 17/09/2019) 
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
%title(hplot,myTitle222)
title(hplot,'log(Kn) std')
grid on

%hplot = subplot(2,2,3);
hplot = subplot(2,3,3);
%himg = imagesc(ESDave,[ESDaveMin,ESDaveMax]);
%imagesc(YAVEssp,[TOPTaveMin,TOPTaveMax]); % Case where ESDave is replaced by optimal temperature (Le Gland, 02/09/2019)
imagesc(YAVEssp,[22,26]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
% Do not cut first and last row in half (Le Gland, 17/09/2019) 
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle223)
grid on

%hplot = subplot(2,2,4);
hplot = subplot(2,3,6);
%%himg = imagesc(ESDstd);
%himg = imagesc(ESDcv,[ESDstdMin,ESDstdMax]);
imagesc(YSTDssp,[TOPTstdMin,TOPTstdMax]); % Case where ESDave is replaced by optimal temperature (Le Gland, 02/09/2019)
%imagesc(YSTDssp,[0.9,1.8]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
% Do not cut first and last row in half (Le Gland, 17/09/2019) 
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle224)
grid on

hplot = subplot(2,3,4);
%imagesc(XYCORssp,[-CorrelationAbsMax +CorrelationAbsMax]);
imagesc(XYCORssp,[-0.6 +0.6]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
% Do not cut first and last row in half (Le Gland, 17/09/2019) 
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle225)
grid on

hplot = subplot(2,3,1);
%imagesc(XYENTssp,[EntropyMin EntropyMax]);
imagesc(PHYssp,[0 max(PHYssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
% Do not cut first and last row in half (Le Gland, 17/09/2019) 
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle226)
grid on
return
