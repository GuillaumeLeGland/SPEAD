function [hfig] = SPEAD_gaussecomodel1D_imagescmodvsobs(PP_obs,CHL_obs,DIN_obs,PON_obs,PPsspcont,CHLsspcont,DINsspcont,PONsspcont,...
    myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages)

%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 

ndepths_obs = size(CHL_obs,1); % ndepths of data
ndepths_mod = size(CHLsspcont,1); % ndepths of model
ndays = size(CHLsspcont,2); % ndays of model

PPmax  = max([PP_obs(:);PPsspcont(:)]);
CHLmax = max([CHL_obs(:);CHLsspcont(:)]);
DINmax = max([DIN_obs(:);DINsspcont(:)]);
PONmax = max([PON_obs(:);PONsspcont(:)]);

%FIGURES:
hfig = figure(fignum);
%..........................................................................
hplot = subplot(4,2,1);
himg = imagesc(PPsspcont,[0 PPmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays],'Xtick',[ndays/12:ndays/12:ndays])
set(hplot,'Ylim',[0.5 ndepths_mod + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'a) Model primary production (mgC/m3/day)')
grid on
hplot = subplot(4,2,3);
himg = imagesc(CHLsspcont,[0 CHLmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays],'Xtick',[ndays/12:ndays/12:ndays])
set(hplot,'Ylim',[0.5 ndepths_mod + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'c) Model chlorophyll concentration (mg/m3)')
grid on
hplot = subplot(4,2,5);
himg = imagesc(DINsspcont,[0 DINmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays],'Xtick',[ndays/12:ndays/12:ndays])
set(hplot,'Ylim',[0.5 ndepths_mod + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'e) Model dissolved inorganic nitrogen (mmol/m3)')
grid on
hplot = subplot(4,2,7);
himg = imagesc(PONsspcont,[0 PONmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays],'Xtick',[ndays/12:ndays/12:ndays])
set(hplot,'Ylim',[0.5 ndepths_mod + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'g) Model particulate organic nitrogen (mmol/m3)')
grid on
%..........................................................................
hplot = subplot(4,2,2);
himg = imagesc(PP_obs,[0 PPmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[0.5 12.5  ],'Xtick',1:12)
set(hplot,'Ylim',[0.5 ndepths_obs+0.5],'Ytick',0.5:50:200.5,'YtickLabel',num2str([0:50:200]'))
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [month]')
title(hplot,'b) Observed primary production (mgC/m3/day)')
grid on
hplot = subplot(4,2,4);
himg = imagesc(CHL_obs,[0 CHLmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[0.5 12.5  ],'Xtick',1:12)
set(hplot,'Ylim',[0.5 ndepths_obs+0.5],'Ytick',0.5:50:200.5,'YtickLabel',num2str([0:50:200]'))
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [month]')
title(hplot,'d) Observed chlorophyll concentration (mg/m3)')
grid on
hplot = subplot(4,2,6);
himg = imagesc(DIN_obs,[0 DINmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[0.5 12.5  ],'Xtick',1:12)
set(hplot,'Ylim',[0.5 ndepths_obs+0.5],'Ytick',0.5:50:200.5,'YtickLabel',num2str([0:50:200]'))
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [month]')
title(hplot,'f) Observed dissolved inorganic nitrogen (mmol/m3)')
grid on
hplot = subplot(4,2,8);
% $$$ plot(tode,ZOO,'-b.')
% $$$ set(hplot,'Ylim',[0 6])
himg = imagesc(PON_obs,[0 PONmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[0.5 12.5  ],'Xtick',1:12)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths_obs+0.5],'Ytick',0.5:50:200.5,'YtickLabel',num2str([0:50:200]'))
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [month]')
title(hplot,'h) Observed particulate organic nitrogen (mmol/m3)')
grid on

