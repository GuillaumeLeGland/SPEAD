function [hfig] = SPEAD_1D_imagescNPZD(NTOTssp,CHLssp,PPssp,ZOOssp,DINssp,PONssp,ndepths,ndays,CHLmin,CHLmax,ZOOmin,ZOOmax,...
    DINmin,DINmax,PONmin,PONmax,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages)

%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 

%FIGURES:
hfig = figure(fignum);
hplot = subplot(2,3,1);
himg = imagesc(NTOTssp,[0 max(NTOTssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'NTOT')
grid on
%..........................................................................
hplot = subplot(2,3,2);
%himg = imagesc(PHYssp,[PHYmin PHYmax]);
himg = imagesc(CHLssp,[CHLmin CHLmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'PHY')
grid on
%..........................................................................
hplot = subplot(2,3,3);
%PP = PHYssp.*MUPssp;
PP = PPssp;
himg = imagesc(PP,[0 max(PP(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'Primary Production')
grid on
%..........................................................................
hplot = subplot(2,3,4);
himg = imagesc(ZOOssp,[ZOOmin ZOOmax]);
%himg = imagesc(ZOOssp,[0 max(ZOOssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'ZOO')
grid on
%..........................................................................
hplot = subplot(2,3,5);
himg = imagesc(DINssp,[DINmin DINmax]);
%himg = imagesc(DINssp,[0 max(DINssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'DIN')
grid on
%..........................................................................
hplot = subplot(2,3,6);
himg = imagesc(PONssp,[PONmin PONmax]);
%himg = imagesc(PONssp,[0 max(PONssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'PON')
grid on
