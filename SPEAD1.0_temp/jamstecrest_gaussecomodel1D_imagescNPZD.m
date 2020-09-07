% Plot article figure with Ntot, phy, PP (PHY*MUP), zoo, DIN and PON (Le
% Gland, (18/11/2019)
%function [hfig] = jamstecrest_gaussecomodel1D_NPZD(NTOTssp,PHYssp,MUPssp,ZOOssp,DINssp,PONssp,fignum,mypackages)
function [hfig] = jamstecrest_gaussecomodel1D_NPZD(NTOTssp,CHLssp,PPssp,ZOOssp,DINssp,PONssp,fignum,mypackages)
global myXtickMarks myXtickLabel myXaxisLabel
global myYtickMarks myYtickLabel myYaxisLabel
global zdepths ndepths deltaz
global t0 deltat ndays nyear tmax 
global PHYmax ZOOmax DINmax PONmax 
global PHYmin ZOOmin DINmin PONmin 

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
hplot = subplot(2,3,2);
% $$$ plot(tode,PHY,'-g.')
% $$$ set(hplot,'Ylim',[0 6])
%himg = imagesc(PHYssp,[PHYmin PHYmax]);
himg = imagesc(CHLssp,[0 max(CHLssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'PHY')
grid on
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
hplot = subplot(2,3,4);
% $$$ plot(tode,ZOO,'-b.')
% $$$ set(hplot,'Ylim',[0 6])
%himg = imagesc(ZOOssp,[ZOOmin ZOOmax]);
himg = imagesc(ZOOssp,[0 max(ZOOssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'ZOO')
grid on
hplot = subplot(2,3,5);
% $$$ plot(tode,DIN,'-r.')
% $$$ set(hplot,'Ylim',[0 6])
% himg = imagesc(DINssp,[DINmin DINmax]);
himg = imagesc(DINssp,[0 max(DINssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'DIN')
grid on
hplot = subplot(2,3,6);
% $$$ plot(tode,PON,'-k.')
% $$$ %%plot(tode,NTOT,'-y.')
% $$$ set(hplot,'Ylim',[0 6])
%himg = imagesc(PONssp,[PONmin PONmax]);
himg = imagesc(PONssp,[0 max(PONssp(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
%%title(hplot,'Ntot')
title(hplot,'PON')
grid on
