function [gca000,hcbar000] = myimagesc2Dhovmoller(AdataHov,myXlims,myYlims,myXtick,myYtick,myXtickLabel,myYtickLabel,cmap,hcmaplims)

%..................................................................
myFontsize = [08]; %For 3x3 subplots.
% $$$ myFontsize = [10]; %For 2x1 subplots.
%..................................................................
myXtickMedio = (myXtick(1:end-1) + myXtick(2:end))/2;
%..................................................................
myXtickAxes001 = myXtickMedio;
myXtickLabelAxes001 = myXtickLabel;
%..................................................................
myXtickAxes002 = myXtick;
%..................................................................
himg001 = imagesc(AdataHov,hcmaplims);
gca001 = gca; %Must be *bofore* colorbar.
gca001position = get(gca001, 'Position');
[hcbar001] = mycolorbar('vertic');
set(hcbar001,'Visible','off');
set(gca001,'Xlim',myXlims,'XTick',myXtickAxes001,'XTickLabel',myXtickLabelAxes001,'Fontsize',myFontsize)
set(gca001,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',[])
%..................................................................
gca002 = axes('Position',gca001position,'Visible','off');
himg002 = imagesc(AdataHov,hcmaplims);
[hcbar002] = mycolorbar('vertic');
set(hcbar002,'Visible','on');
set(gca002,'Xlim',myXlims,'Xtick',myXtickAxes002,'XtickLabel',[])
set(gca002,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel,'Fontsize',myFontsize)
%..................................................................
% $$$ set(get(gca001,'XLabel'),'String','myXlabel')
% $$$ set(get(gca002,'YLabel'),'String','myXlabel')
%..................................................................
colormap(cmap)
grid on
%..................................................................
gca000(1) = gca001;
gca000(2) = gca002;
hcbar000 = hcbar002;
%..................................................................
return


colormap(cmap)
set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',[])
set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)

set(hcbar002,'Ylim',hcmaplims,'Ytick',hcmapYtick001,'YtickLabel',hcmapYtickLabel001)
% $$$ xlabel('Time')
ylabel('Latitude')
htitle001 = title('a) Primary Production')
grid on


