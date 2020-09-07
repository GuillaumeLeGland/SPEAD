function [hfig] = jamstecrest_gaussecomodel1D_imagescuptakerates(MUP,MUZ,Array2D,NTOT,fignum,mypackages)
%...................................................................................
global myTitle001 myTitle002 myTitle003 myTitle004 
global myXtickMarks myXtickLabel myXaxisLabel
global myYtickMarks myYtickLabel myYaxisLabel
global zdepths ndepths deltaz
global t0 deltat ndays nyear tmax 
%...................................................................................
global MUPmin MUPmax
global MUZmin MUZmax 
%...................................................................................
global minArray2D maxArray2D 
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
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'Uptake rate [d-1]')
grid on
hplot = subplot(2,2,2);
himg = imagesc(MUZ,[MUZmin MUZmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'Grazing rate [d-1]')
grid on
hplot = subplot(2,2,3);
himg = imagesc(Array2D,[minArray2D-sqrt(eps) maxArray2D]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,myTitle003)
grid on
hplot = subplot(2,2,4);
himg = imagesc(NTOT,[0 max(NTOT(:))]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Xlim',[1 ndays  ],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
%set(hplot,'Ylim',[1 ndepths],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
set(hplot,'Ylim',[0.5 ndepths + 0.5],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
ylabel(hplot,myYaxisLabel)
xlabel(hplot,'time [days]')
title(hplot,'NTOT')
grid on
return
