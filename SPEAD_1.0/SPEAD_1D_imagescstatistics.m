function [hfig] = SPEAD_1D_imagescstatistics(XAVEssp,XSTDssp,YAVEssp,YSTDssp,XYCORssp,PHYssp,ndepths,ndays,myXtickMarks,myXtickLabel,...
    myYtickMarks,myYtickLabel,myYaxisLabel,logESDaveMax,logESDaveMin,logESDstdMax,logESDstdMin,TOPTaveMax,TOPTaveMin,TOPTstdMax,TOPTstdMin,...
    CorrelationAbsMax,fignum,mypackages)
%...................................................................................
%global myXtickMarks myXtickLabel
%global myYtickMarks myYtickLabel myYaxisLabel
%global ndepths ndays 
%...................................................................................
%global logESDaveMax logESDstdMax 
%global logESDaveMin logESDstdMin  
%global   TOPTaveMax   TOPTstdMax
%global   TOPTaveMin   TOPTstdMin
%global         CorrelationAbsMax 
%...................................................................................
myTitle221 = 'log(ESD) ave';
myTitle222 = 'log(ESD) std';
myTitle223 = 'Topt ave';
myTitle224 = 'Topt std';
myTitle225 = 'Correlation';
myTitle226 = 'Total biomass';

%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 

%FIGURES:
hfig = figure(fignum);
hplot = subplot(2,3,1);
imagesc(PHYssp,[0 max(PHYssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel) 
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle226)
grid on
%..........................................................................
hplot = subplot(2,3,2);
%himg = imagesc(XAVEssp,[min(XAVEssp(:)),max(XAVEssp(:))]);
himg = imagesc(XAVEssp,[logESDaveMin,logESDaveMax]);
%himg = imagesc(XAVEssp,[-0.8,+0.6]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel) 
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle221)
grid on
%..........................................................................
hplot = subplot(2,3,3);
%himg = imagesc(YAVEssp,[min(YAVEssp(:)),max(YAVEssp(:))]);
himg = imagesc(YAVEssp,[TOPTaveMin,TOPTaveMax]);
%imagesc(YAVEssp,[22,26]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle223)
grid on
%..........................................................................
hplot = subplot(2,3,4);
%imagesc(XYCORssp,[-max(abs(XYCORssp(:))) max(abs(XYCORssp(:)))]);
imagesc(XYCORssp,[-CorrelationAbsMax +CorrelationAbsMax]);
%imagesc(XYCORssp,[-0.6 +0.6]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle225)
grid on
%..........................................................................
hplot = subplot(2,3,5);
himg = imagesc(XSTDssp,[logESDstdMin,logESDstdMax]);
%himg = imagesc(XSTDssp,[min(XSTDssp(:)),max(XSTDssp(:))]);
%himg = imagesc(XSTDssp,[0.3,0.6]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle222)
grid on
%..........................................................................
hplot = subplot(2,3,6);
himg = imagesc(YSTDssp,[TOPTstdMin,TOPTstdMax]);
%himg = imagesc(YSTDssp,[min(YSTDssp(:)),max(YSTDssp(:))]);
%imagesc(YSTDssp,[0.9,1.8]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle224)
grid on
%..........................................................................
return
