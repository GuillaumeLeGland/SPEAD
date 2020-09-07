function [ ] = jamstecrest_imagescforcings(temp,parz0,par2D,mld,kz,fignum,mypackages)
% Plots the forcing of the model (Temperature, irrradiance, mixed layer
% depth and vertical diffusivity) (Le Gland, 29/10/2019)

%MY PACKAGES FOR PLOTING:
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales; 

hfig = figure(fignum);

% Temperature subplot
hplot = subplot(2,2,1);
tmin = floor(min(temp(:)));
tmax = ceil(max(temp(:)));
himg = imagesc(temp,[tmin,tmax]);
hcbar = colorbar_funhan(verticales);
set(hplot,'Ylim',[0.5 20.5],'Ytick',[0.5:5:20.5],'YtickLabel',num2str([0:50:200]'))
%set(hplot,'Xlim',[1 360  ],'Xtick',[15:30:345],'XtickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'])
set(hplot,'Xlim',[0 360  ],'Xtick',[0:30:360],'XtickLabel',num2str([0:30:360]'))
ylabel(hplot,'Depth [m]')
xlabel(hplot,'Time [days]')
grid on
title(['a) Temperature [' char(176) 'C]'])

% PAR subplot in 2D
hplot = subplot(2,2,2);
plot(parz0,'k-','Linewidth',3)
set(hplot,'Ylim',[0 150.0],'Ytick',[0:25:150],'YtickLabel',num2str([0:25:150]'))
set(hplot,'Xlim',[0 360  ],'Xtick',[0:30:360],'XtickLabel',num2str([0:30:360]'))
ylabel(hplot,'PAR [W.m^{-2}]')
xlabel(hplot,'Time [days]')
grid on
title('b) Photosynthetically Available Radiation at surface (W.m^{-2})')

% PAR subplot in 2D
hplot = subplot(2,2,3);
logpar = log10(par2D);
logparmin = 0.1*floor(10*min(logpar(:)));
logparmax = 0.1*ceil(10*max(logpar(:)));
himg = imagesc(logpar,[logparmin,logparmax]);
%himg = imagesc(par);
hcbar = colorbar_funhan(verticales);
%hcbar = colorbar;
caxis([-1 logparmax])
L = [0.03 0.1 0.3 1 3 10 30 100];
l = (log10(L)+1)-1;
set(hcbar,'Ytick',l,'YTicklabel',L);
set(hplot,'Ylim',[0.5 20.5],'Ytick',[0.5:5:20.5],'YtickLabel',num2str([0:50:200]'))
set(hplot,'Xlim',[0 360  ],'Xtick',[0:30:360],'XtickLabel',num2str([0:30:360]'))
ylabel(hplot,'Depth [m]')
xlabel(hplot,'Time [days]')
grid on
title('c) Photosynthetically Available Radiation at depth (W.m^{-2})')

% Vertical diffusivity(kz) subplot
hplot = subplot(2,2,4);
logkz = log10(kz);
logkzmin = floor(min(logkz(:)));
logkzmax = ceil(max(logkz(:)));
%kzmin = 10^(floor(min(log10(kz(:)))));
%kzmax = 10^(ceil(max(log10(kz(:)))));
himg = imagesc(logkz,[logkzmin,logkzmax]);
hcbar = colorbar_funhan(verticales);
%hcbar = colorbar;
L = [0.1 1 10 100 1000 10000];
l = (log10(L)+1)-1;
set(hcbar,'Ytick',l,'YTicklabel',L);
% 
set(hplot,'Ylim',[1 21],'Ytick',[1:5:21],'YtickLabel',num2str([0:50:200]'))
%set(hplot,'Xlim',[1 360  ],'Xtick',[15:30:345],'XtickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'])
set(hplot,'Xlim',[0 360  ],'Xtick',[0:30:360],'XtickLabel',num2str([0:30:360]'))
ylabel(hplot,'Depth [m]')
xlabel(hplot,'Time [days]')
hold on
plot((mld+10)/10,'k-','LineWidth',3)
grid on
title('d) Vertical diffusivity K_z (m^2.d^{-1}) and mixed layer depth (m)')


end

