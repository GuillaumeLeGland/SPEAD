function [hfig] = jamstecrest_loopruns_imagescsubplots_4x4(aveMUP,aveMUZ,aveFPHYT,aveFZOO,avePHYT,aveZOO,aveDIN,avePON,aveESDphyAve,aveEindex,fignum,mypackages)
%....................................................................
global myYtickMarks myYtickLabel myYaxisLabel 
global myXtickMarks myXtickLabel myXaxisLabel 
%....................................................................

%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 

%FIGURES:
%====================================================================
%....................................................................
[msize,nsize] = size(aveMUP);
%....................................................................
%====================================================================
%%hfig = figure('Visible','off','PaperPosition',[0 0 6 4],'PaperSize',[6 4]);
hfig = figure(fignum);
%....................................................................
subplot(3,3,1)
imagesc(aveMUP)
set(gca,'Xlim',[1 msize],'Xtick',myXtickMarks,'XtickLabel','')
set(gca,'Ylim',[1 nsize],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
%%xlabel(myXaxisLabel)
ylabel(myYaxisLabel)
colorbar_funhan(verticales)
%%title('muphy')
title('Productivity Phy')
axis square 
grid on
subplot(3,3,2)
imagesc(aveFPHYT)
set(gca,'Xlim',[1 msize],'Xtick',myXtickMarks,'XtickLabel','')
set(gca,'Ylim',[1 nsize],'Ytick',myYtickMarks,'YtickLabel','')
%%xlabel(myXaxisLabel)
%%ylabel(myYaxisLabel)
colorbar_funhan(verticales)
%%title('Fphy')
title('Production Phy')
axis square 
grid on
subplot(3,3,3)
imagesc(avePHYT)
set(gca,'Xlim',[1 msize],'Xtick',myXtickMarks,'XtickLabel','')
set(gca,'Ylim',[1 nsize],'Ytick',myYtickMarks,'YtickLabel','')
%%xlabel(myXaxisLabel)
%%ylabel(myYaxisLabel)
colorbar_funhan(verticales)
%%title('phy')
title('Biomass Phy')
axis square 
grid on
%....................................................................
subplot(3,3,1+3)
imagesc(aveMUZ)
set(gca,'Xlim',[1 msize],'Xtick',myXtickMarks,'XtickLabel','')
set(gca,'Ylim',[1 nsize],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
%%xlabel(myXaxisLabel)
ylabel(myYaxisLabel)
colorbar_funhan(verticales)
%%title('muzoo')
title('Productivity Zoo')
axis square 
grid on
subplot(3,3,2+3)
imagesc(aveFZOO)
set(gca,'Xlim',[1 msize],'Xtick',myXtickMarks,'XtickLabel','')
set(gca,'Ylim',[1 nsize],'Ytick',myYtickMarks,'YtickLabel','')
%%xlabel(myXaxisLabel)
%%ylabel(myYaxisLabel)
colorbar_funhan(verticales)
%%title('Fzoo')
title('Production Zoo')
axis square 
grid on
subplot(3,3,3+3)
imagesc(aveZOO)
set(gca,'Xlim',[1 msize],'Xtick',myXtickMarks,'XtickLabel','')
set(gca,'Ylim',[1 nsize],'Ytick',myYtickMarks,'YtickLabel','')
%%xlabel(myXaxisLabel)
%%ylabel(myYaxisLabel)
colorbar_funhan(verticales)
%%title('zoo')
title('Biomass Zoo')
axis square 
grid on
%....................................................................
subplot(3,3,1+6)
imagesc(aveEindex)
set(gca,'Xlim',[1 msize],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(gca,'Ylim',[1 nsize],'Ytick',myYtickMarks,'YtickLabel',myYtickLabel)
xlabel(myXaxisLabel)
ylabel(myYaxisLabel)
colorbar_funhan(verticales)
title('Diversity')
%%title('exp(Shannon)')
axis square 
grid on
subplot(3,3,2+6)
%%imagesc(avePON)
imagesc(aveESDphyAve)
set(gca,'Xlim',[1 msize],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(gca,'Ylim',[1 nsize],'Ytick',myYtickMarks,'YtickLabel','')
xlabel(myXaxisLabel)
%%ylabel(myYaxisLabel)
colorbar_funhan(verticales)
%%title('PON')
title('ESD ave')
axis square 
grid on
subplot(3,3,3+6)
imagesc(aveDIN)
set(gca,'Xlim',[1 msize],'Xtick',myXtickMarks,'XtickLabel',myXtickLabel)
set(gca,'Ylim',[1 nsize],'Ytick',myYtickMarks,'YtickLabel','')
xlabel(myXaxisLabel)
%%ylabel(myYaxisLabel)
colorbar_funhan(verticales)
title('DIN')
axis square 
grid on
%....................................................................
%====================================================================
%********************************************************************
return

