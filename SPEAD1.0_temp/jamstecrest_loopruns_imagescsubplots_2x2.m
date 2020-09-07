function [ ] = jamstecrest_loopruns_imagescsubplots_2x2(A001,A002,A003,A004,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages)

%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 

%FIGURES:
%====================================================================
figure(fignum)
%....................................................................
subplot(2,2,1)
imagesc(A001)
colorbar_funhan(verticales)
title(myTitle001)
%....................................................................
subplot(2,2,2)
imagesc(A002)
colorbar_funhan(verticales)
title(myTitle002)
%....................................................................
subplot(2,2,3)
imagesc(A003)
colorbar_funhan(verticales)
title(myTitle003)
%....................................................................
subplot(2,2,4)
imagesc(A004)
colorbar_funhan(verticales)
title(myTitle004)
%....................................................................
%====================================================================
%********************************************************************
return

