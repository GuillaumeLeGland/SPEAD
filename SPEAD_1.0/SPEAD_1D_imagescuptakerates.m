function [hfig] = SPEAD_1D_imagescuptakerates(MUP,MUZ,NTOT,ndepths,ndays,MUPmin,MUPmax,MUZmin,MUZmax,...
    myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages)
%...................................................................................

%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 

%FIGURES:
hfig = figure(fignum);
hplot = subplot(2,2,1);
himg = imagesc(MUP,[MUPmin MUPmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'Uptake rate [d-1]')
grid on
%..........................................................................
hplot = subplot(2,2,2);
himg = imagesc(MUZ,[MUZmin MUZmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'Grazing rate [d-1]')
grid on
%..........................................................................
hplot = subplot(2,2,4);
himg = imagesc(NTOT,[0 max(NTOT(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'NTOT')
grid on
return
