function [ax2,himg2,hbar2] = jamstecrest_gaussecomodel1D_interpimagesc(A,Amin,Amax,ndays,SizeFontTicks,myTitle,mypackages)
global myYlims1m myYtick1m myYtickLabel1m 
global myXlims1m myXtick1m myXtickLabel1m
global monthlimits

%===================================================================================
%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 
%===================================================================================
%FIGURES:
%...................................................................................
% $$$ Amin = min(A(:));
% $$$ Amax = max(A(:));
%...................................................................................
Adel = (Amax-Amin)/4;
%...................................................................................
[AI] = jamstecrest_interp2D(A,ndays);
%...................................................................................
himg1 = imagesc(AI,[Amin Amax]);
ax1 = gca; %Must be before the "colorbar".
ax1position = get(ax1, 'Position');
hbar1 = colorbar_funhan(horizontal);
set(ax1,'Ylim',myYlims1m,'Ytick',myYtick1m,'YTickLabel',myYtickLabel1m,'FontSize',SizeFontTicks)
set(ax1,'Xlim',myXlims1m,'XTick',myXtick1m,'XTickLabel',myXtickLabel1m,'Fontsize',SizeFontTicks)
set(hbar1,'Visible','off'); %when using ax2 
%...................................................................................
ax2 = axes('Position', ax1position,'Visible', 'off');
himg2 = imagesc(AI,[Amin Amax]);
hbar2 = colorbar_funhan(horizontal);
set(ax2,'Ylim',myYlims1m,'Ytick',myYtick1m,'YTickLabel',[])
set(ax2,'Xlim',myXlims1m,'XTick',monthlimits,'XTickLabel',[])
%%set(hbar2,'Xtick',[Amin:Adel:Amax],'XtickLabel',[Amin:Adel:Amax],'Fontsize',SizeFontTicks)
set(hbar2,'Fontsize',SizeFontTicks)
title(myTitle)
grid on
%...................................................................................
%===================================================================================
return
