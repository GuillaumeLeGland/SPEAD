% Trait distribution at each season and different depths, in the continuous
% and discrete cases (Le Gland, 18/02/2020)
function [] = SPEAD_1D_distribution(PHYTsspcont,logESDphysspAveCont,logESDphysspStdCont,TOPTphysspAveCont,TOPTphysspStdCont,physspCorCont,PHYsspdisc3D,xrng,yrng,itemp,DINsspdisc,fignum,mypackages) 

%============================================================================
%............................................................................
Fontsize = 20;
Markersize = 5;
%............................................................................
xIndex = 1:1:length(xrng);
yIndex = 1:1:length(yrng);
%............................................................................
ymin = 0;
%ymax = 0.16;
%............................................................................

%============================================================================

% colorbar_funhan = mypackages.colorbar;
% verticales = mypackages.verticales;

%ymax = max(PHYsspdisc3D(:)); % Maximum concentration of a single phenotype
%ymin = (1/128)*ymax;
%ymax = 0.01;
ymax = max(PHYsspdisc3D(:));

% Recompute nxphy and nyphy to adapt to the size of the discrete model (Le Gland, 20/07/2020)
nxphy = length(xrng);
nyphy = length(yrng);

% Grid for the ellipse plots
[xmesh, ymesh] = meshgrid(1:0.1:25, 1:0.1:25);
xval = -2.5 +  4*(xmesh-1)/(nxphy-1); % Values of trait x are between -2.5 and +1.5
yval = 18.0 + 12*(ymesh-1)/(nyphy-1); % Values of trait y are between 18.0 and 30.0

% Depth and day of plot
deparr = 1:5:11;
dayarr =71:90:341;

% Color map
myColorMap = [1 1 1; 0 0 0.5; 0 0 0.5; 0 0.85 0.85; 0 0.85 0.85; 0 1 0; 0.9 0.9 0; 1 0.5 0; 0.9 0 0; 0.4 0 0];

Titles = ['a) Winter, surface'; 'b) Spring, surface'; 'c) Summer, surface'; 'd) Autumn, surface'; ...
          'e) Winter, 50 m   '; 'f) Spring, 50 m   '; 'g) Summer, 50 m   '; 'h) Autumn, 50 m   '; ...
          'i) Winter, 100 m  '; 'j) Spring, 100 m  '; 'k) Summer, 100 m  '; 'h) Autumn, 100 m  '];

figure(fignum)
      
for jday = 1:4
    for jk = 1:3
        
        dep = deparr(jk);
        day = dayarr(jday);
        nplot = jday + 4*(jk-1);

        hplot = subplot(3,4,nplot);
        % Color plot of the discrete model
        % White from 0 to exp(-4.5). Correspond to 3 std from a maximum of 1 (Le Gland, 20/07/2020)
        imagesc(squeeze(log(PHYsspdisc3D(dep,day,:,:)/ymax))',[-5.0, 0])
        
        colormap(myColorMap);
        
        xtick = [1:3:nxphy];
        set(hplot,'Xlim',[0.5 length(xrng) + 0.5],'Xtick',xtick,'XtickLabel',0.01*round(100*exp([-2.5:12/(nxphy-1):1.5])))
        set(hplot,'Ylim',[0.5 length(yrng) + 0.5],'Ytick',[1:(nyphy-1)/6:nyphy],'YtickLabel',[18; 20; 22; 24; 26; 28; 30])
        % set(hplot,'Xlim',[0.5 length(xrng) + 0.5],'Xtick',xtick,'XtickLabel',0.01*round(100*exp([-1:6/(nxphy-1):1])))
        % set(hplot,'Ylim',[0.5 length(yrng) + 0.5],'Ytick',[1:(nyphy-1)/6:nyphy],'YtickLabel',[21; 22; 23; 24; 25; 26; 27])
        
        set(gca,'YDir','normal')
        xlabel(hplot,'Half-saturation [mmol.m^{-3}]')
        ylabel(hplot,'Optimal temperature [°C]')
        title(hplot,Titles(nplot,:))
        hold on
        % Ellipse plot of the continuous model
        R = physspCorCont(dep,day);
        xm = logESDphysspAveCont(dep,day);
        ym = TOPTphysspAveCont(dep,day);
        xstd = logESDphysspStdCont(dep,day);
        ystd = TOPTphysspStdCont(dep,day);
        ellipse = -1/(2*(1-R^2)) * ( ((xval-xm)/xstd).^2 - 2*R*(xval-xm).*(yval-ym)/(xstd*ystd) + ((yval-ym)/ystd).^2 );
        % factor to express ellipses relatively to the all-time and all-depth maximum
        aggmax = PHYTsspcont(dep,day) * (xrng(2)-xrng(1)) * (yrng(2)-yrng(1)) * (1/(2*pi)) * 1/(xstd*ystd*sqrt(1-R^2));
        fac = ymax/aggmax;
        
        % Contours should represent 1 std, 2 std and 3 std away from the peak of
        % the distribution. These correspond to 60.7\%, 13.5\% and 1.1\% of the maximum.
        contour(xmesh,ymesh,ellipse + 0.5 - log(fac),[0.0 0.0],'k-','LineWidth',5)
        contour(xmesh,ymesh,ellipse + 2.0 - log(fac),[0.0 0.0],'k-','LineWidth',3)
        contour(xmesh,ymesh,ellipse + 4.5 - log(fac),[0.0 0.0],'k-','LineWidth',2)
        
        plot([1+((nxphy-1)/4)*(2.5+log(DINsspdisc(dep,day))),1+((nxphy-1)/4)*(2.5+log(DINsspdisc(dep,day)))],[1 nxphy],'k--')
        % Plot the expected optimal temperature (T+2) instead of the real environment temperature (T)
        plot([1 nyphy],[1+((nyphy-1)/12)*(itemp(dep,day)-16) 1+((nyphy-1)/12)*(itemp(dep,day)-16)],'k--')
        hold off
        
    end
end


return