function [hfig] = jamstecrest_gaussecomodel1D_imagescbiomasses(PHYssp,ZOOssp,DINssp,PONssp,fignum,mypackages)
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
hplot = subplot(2,2,1);
% $$$ plot(tode,PHY,'-g.')
% $$$ set(hplot,'Ylim',[0 6])
himg = imagesc(PHYssp,[PHYmin PHYmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'PHY')
grid on
hplot = subplot(2,2,2);
% $$$ plot(tode,ZOO,'-b.')
% $$$ set(hplot,'Ylim',[0 6])
himg = imagesc(ZOOssp,[ZOOmin ZOOmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'ZOO')
grid on
hplot = subplot(2,2,3);
% $$$ plot(tode,DIN,'-r.')
% $$$ set(hplot,'Ylim',[0 6])
himg = imagesc(DINssp,[DINmin DINmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'DIN')
grid on
hplot = subplot(2,2,4);
% $$$ plot(tode,PON,'-k.')
% $$$ %%plot(tode,NTOT,'-y.')
% $$$ set(hplot,'Ylim',[0 6])
himg = imagesc(PONssp,[PONmin PONmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
%%title(hplot,'Ntot')
title(hplot,'PON')
grid on
