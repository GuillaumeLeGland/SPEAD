function [hcbar] = SPEAD_1D_subplotesting(A,Amin,Amax,fignum,mypackages)

%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 

%FIGURES:
[msize,nsize] = size(A);

figure(fignum + 1);
hplot = subplot(2,2,1);
himg = imagesc(A,[Amin Amax]);
hcbar = colorbar(verticales);
ylabel(hplot,'Y axis')
xlabel(hplot,'X axis')
set(hplot,'Xlim',[0 nsize],'Xtick',[0:64:nsize])
set(hplot,'Ylim',[0 msize],'Ytick',[0:32:nsize])
title(hplot,'PLOT 1')
hplot = subplot(2,2,2);
himg = imagesc(A,[Amin Amax]);
hcbar = colorbar(horizontal);
ylabel(hplot,'Y axis')
xlabel(hplot,'X axis')
set(hplot,'Xlim',[0 nsize],'Xtick',[0:64:nsize])
set(hplot,'Ylim',[0 msize],'Ytick',[0:32:nsize])
title(hplot,'PLOT 2')

figure(fignum + 2);
hplot = subplot(2,2,1);
himg = imagesc(A,[Amin Amax]);
hcbar = colorbar_funhan(verticales);
ylabel(hplot,'Y axis')
xlabel(hplot,'X axis')
set(hplot,'Xlim',[0 nsize],'Xtick',[0:64:nsize])
set(hplot,'Ylim',[0 msize],'Ytick',[0:32:nsize])
title(hplot,'PLOT 1')
hplot = subplot(2,2,2);
himg = imagesc(A,[Amin Amax]);
hcbar = colorbar_funhan(horizontal);
ylabel(hplot,'Y axis')
xlabel(hplot,'X axis')
set(hplot,'Xlim',[0 nsize],'Xtick',[0:64:nsize])
set(hplot,'Ylim',[0 msize],'Ytick',[0:32:nsize])
title(hplot,'PLOT 2')

figure(fignum + 3);
hplot = subplot_funhan(2,2,1);
himg = imagesc(A,[Amin Amax]);
hcbar = colorbar_funhan(verticales); 
ylabel(hplot,'Y axis')
xlabel(hplot,'X axis')
set(hplot,'Xlim',[0 nsize],'Xtick',[0:64:nsize])
set(hplot,'Ylim',[0 msize],'Ytick',[0:32:nsize])
title(hplot,'PLOT 1')
hplot = subplot_funhan(2,2,2);
himg = imagesc(A,[Amin Amax]);
hcbar = colorbar_funhan(horizontal);
ylabel(hplot,'Y axis')
xlabel(hplot,'X axis')
set(hplot,'Xlim',[0 nsize],'Xtick',[0:64:nsize])
set(hplot,'Ylim',[0 msize],'Ytick',[0:32:nsize])
title(hplot,'PLOT 2')

return

