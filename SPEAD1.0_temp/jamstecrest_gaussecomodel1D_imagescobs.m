% Plot article figure with Ntot, phy, PP (PHY*MUP), zoo, DIN and PON (Le
% Gland, (18/11/2019)
function [hfig] = jamstecrest_gaussecomodel1D_imagescobs(CHL_obs,PP_obs,DIN_obs,PON_obs,fignum,mypackages)
global myXtickMarks myXtickLabel myXaxisLabel
global myYtickMarks myYtickLabel myYaxisLabel
%global zdepths ndepths deltaz
%global t0 deltat ndays nyear tmax 
%global PHYmax ZOOmax DINmax PONmax 
%global PHYmin ZOOmin DINmin PONmin 

%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 

ndepths = size(CHL_obs,1); % ndepths of data

%FIGURES:
hfig = figure(fignum);
hplot = subplot(2,2,1);
himg = imagesc(CHL_obs,[0 max(CHL_obs(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[0.5 12.5  ])
set(hplot,'Ylim',[0.5 ndepths + 0.5])
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'Chl (mg/m3)')
grid on
hplot = subplot(2,2,2);
% $$$ plot(tode,PHY,'-g.')
% $$$ set(hplot,'Ylim',[0 6])
himg = imagesc(PP_obs,[0 max(PP_obs(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[0.5 12.5  ])
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5])
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'Primary Production (mmolN/m2/day)')
grid on
hplot = subplot(2,2,3);
himg = imagesc(DIN_obs,[0 max(DIN_obs(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[0.5 12.5  ])
set(hplot,'Ylim',[0.5 ndepths + 0.5])
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'Dissolved Inorganic Nitrogen (mmol/m3)')
grid on
hplot = subplot(2,2,4);
% $$$ plot(tode,ZOO,'-b.')
% $$$ set(hplot,'Ylim',[0 6])
himg = imagesc(PON_obs,[0 max(PON_obs(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[0.5 12.5  ])
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5])
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'Particulate Organic Nitrogen (mmol/m3)')

