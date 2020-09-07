% Trait distribution at each season and different depths, in the continuous
% and discrete cases (Le Gland, 18/02/2020)
function [] = jamstecrest_distribution(PHYTsspcont,logESDphysspAveCont,logESDphysspStdCont,TOPTphysspAveCont,TOPTphysspStdCont,physspCorCont,PHYsspdisc3D,xrng,yrng,itemp,DINsspdisc,fignum,mypackages) 

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
%xval = -1 + 2*(xmesh-1)/(nxphy-1); % Values of trait x are between -1.0 and +1.0
%yval = 21.0 + 6*(ymesh-1)/(nyphy-1); % Values of trait y are between 21.0 and 27.0

% Depth and day of plot
deparr = 1:5:11;
dayarr =71:90:341;

% Color map
%myColorMap = jet(13);
%myColorMap(1,:) = 1;
%myColorMap = [1 1 1; jet(8)];
%myColorMap = [1 1 1; jet(9)];
%myColorMap = [1 1 1; 0 0 0.5; 0 0 1; 0 0.5 1; 0 1 1; 0 1 0; 1 1 0; 1 0.5 0; 1 0 0; 0.5 0 0];
%Reduce the number of different blues
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
        % imagesc(squeeze(log2(PHYsspdisc3D(dep,day,:,:)))',[log2(ymin), log2(ymax)])
        % imagesc(squeeze(log2(PHYsspdisc3D(dep,day,:,:)/ymax))',[-log2(100), 0])
        
        % White from 0 to exp(-4.5). Correspond to 3 std from a maximum of 1 (Le Gland, 20/07/2020)
        imagesc(squeeze(log(PHYsspdisc3D(dep,day,:,:)/ymax))',[-5.0, 0])
        
        colormap(myColorMap);
        % hcbar = colorbar_funhan(verticales);
        %set(hplot,'Xlim',[0.5 length(xrng) + 0.5],'Xtick',[4 10 16 22],'XtickLabel',[exp(-2); exp(-1); exp(0); exp(1)])
        %xtick = [2.629,6.961,11.293,15.625,19.957,24.289];
        %set(hplot,'Xlim',[0.5 length(xrng) + 0.5],'Xtick',xtick,'XtickLabel',[0.125; 0.25; 0.5; 1.0; 2.0; 4.0])
        %set(hplot,'Ylim',[0.5 length(yrng) + 0.5],'Ytick',[1 5 9 13 17 21 25],'YtickLabel',[18; 20; 22; 24; 26; 28; 30])
        
        % Be more "honest" aboout the values of half-saturation for each
        % phenotype of the discrete model (20/07/2020)
        %xtick = [1:(nxphy-1)/8:nxphy];
        xtick = [1:3:nxphy];
        set(hplot,'Xlim',[0.5 length(xrng) + 0.5],'Xtick',xtick,'XtickLabel',0.01*round(100*exp([-2.5:12/(nxphy-1):1.5])))
        set(hplot,'Ylim',[0.5 length(yrng) + 0.5],'Ytick',[1:(nyphy-1)/6:nyphy],'YtickLabel',[18; 20; 22; 24; 26; 28; 30])
        % set(hplot,'Xlim',[0.5 length(xrng) + 0.5],'Xtick',xtick,'XtickLabel',0.01*round(100*exp([-1:6/(nxphy-1):1])))
        % set(hplot,'Ylim',[0.5 length(yrng) + 0.5],'Ytick',[1:(nyphy-1)/6:nyphy],'YtickLabel',[21; 22; 23; 24; 25; 26; 27])
        
        set(gca,'YDir','normal')
        xlabel(hplot,'Half-saturation [mmol.m^{-3}]')
        ylabel(hplot,'Optimal temperature [°C]')
        %title(hplot,'a) Winter, surface')
        title(hplot,Titles(nplot,:))
        hold on
        % Ellipse plot of the continuous model
        R = physspCorCont(dep,day);
        xm = logESDphysspAveCont(dep,day);
        ym = TOPTphysspAveCont(dep,day);
        xstd = logESDphysspStdCont(dep,day);
        ystd = TOPTphysspStdCont(dep,day);
        ellipse = -1/(2*(1-R^2)) * ( ((xval-xm)/xstd).^2 - 2*R*(xval-xm).*(yval-ym)/(xstd*ystd) + ((yval-ym)/ystd).^2 );
        %plot(logESDphysspAveCont(1,60),TOPTphysspAveCont(1,60),'k+')
        % factor to express ellipses relatively to the all-time and all-depth maximum
        % fac = ymax/max(PHYsspdisc3D(dep,day,:)); (not the right density factor, should be based on aggregated model itself
        aggmax = PHYTsspcont(dep,day) * (xrng(2)-xrng(1)) * (yrng(2)-yrng(1)) * (1/(2*pi)) * 1/(xstd*ystd*sqrt(1-R^2));
        fac = ymax/aggmax;
        
        %contour(xmesh,ymesh,ellipse + log(5) - log(fac),[0.0 0.0],'k-','LineWidth',3)
        %contour(xmesh,ymesh,ellipse + log(20) - log(fac),[0.0 0.0],'k-','LineWidth',2)
        %contour(xmesh,ymesh,ellipse + log(100)- log(fac),[0.0 0.0],'k-','LineWidth',1)
        
        % Contours should represent 1 std, 2 std and 3 std away from the peak of
        % the distribution. These correspond to 60.7\%, 13.5\% and 1.1\% of the maximum.
        contour(xmesh,ymesh,ellipse + 0.5 - log(fac),[0.0 0.0],'k-','LineWidth',5)
        contour(xmesh,ymesh,ellipse + 2.0 - log(fac),[0.0 0.0],'k-','LineWidth',3)
        contour(xmesh,ymesh,ellipse + 4.5 - log(fac),[0.0 0.0],'k-','LineWidth',2)
        
        plot([1+((nxphy-1)/4)*(2.5+log(DINsspdisc(dep,day))),1+((nxphy-1)/4)*(2.5+log(DINsspdisc(dep,day)))],[1 nxphy],'k--')
        %plot([1 7],[1+0.5*(itemp(dep,day)-18) 1+0.5*(itemp(dep,day)-18)],'k--')
        % Plot the expected optimal temperature (T+2) instead of the real
        % environment temperature (T) (Le Gland, 20/07/2020)
        plot([1 nyphy],[1+((nyphy-1)/12)*(itemp(dep,day)-16) 1+((nyphy-1)/12)*(itemp(dep,day)-16)],'k--')
        hold off
        
        %if jk == 3 && jday == 1
        %    c = colorbar('southoutside');
        %end
        
    end
end


return

%============================================================================
%............................................................................
%Choose locations here (Le Gland, 19/11/2019)
%depth = [6,1,11];
%SamplingDays = [135,225,300];
nsamplingdays = length(SamplingDays);
%cmap = jet(nsamplingdays);
%cmap = ['g';'r';'b'];
cmap = [[0 1 0];[1 0 0];[0 0 1]];
nxphy = length(xrng);
nyphy = length(yrng);
xdel = (xrng(end)-xrng(1))/(nxphy-1);
ydel = (yrng(end)-yrng(1))/(nyphy-1);
phy_cont = ones(nsamplingdays,nxphy,nyphy)*nan;
phy_disc = ones(nsamplingdays,nxphy,nyphy)*nan;
%............................................................................
for j = 1:nsamplingdays 
    phytot = PHYTsspcont(depth(j),SamplingDays(j));
    xave = logESDphysspAveCont(depth(j),SamplingDays(j));
    xsig = logESDphysspStdCont(depth(j),SamplingDays(j));
    yave = TOPTphysspAveCont(depth(j),SamplingDays(j));
    ysig = TOPTphysspStdCont(depth(j),SamplingDays(j));
    xycor = physspCorCont(depth(j),SamplingDays(j));
    xydet = xsig^2 * ysig^2 - (xsig*ysig*xycor)^2;
    xymat = [xsig^2, xsig*ysig*xycor; xsig*ysig*xycor, ysig^2];
    %xyinvmat = inv(xymat);
    PDF = ones(nxphy,nyphy)*nan;
    for jx = 1:nxphy
        for jy = 1:nyphy
            xyvec = [xrng(jx) - xave; yrng(jy) - yave]; % Use xaxis and yaxis instead of creating new variables (Le Gland, 17/09/2019)
            % PDF(jx,jy) = (1.0 / (sqrt(xydet) * 2 * pi)) * exp( - (1/2) * xyvec' * xyinvmat * xyvec);
            % Use xymat \ xyvec instead of xyinvmat * xyvec for efficienncy (Le Gland, 19/11/2019)
            PDF(jx,jy) = (1.0 / (sqrt(xydet) * 2 * pi)) * exp( - (1/2) * xyvec' * (xymat \ xyvec));
            % PDF(jx,jy) = (1.0 / (xydet^(1/4) * 2 * sqrt(pi))) * exp( - (1/2) * xyvec' * (xymat \ xyvec));
        end
    end
    phy_cont(j,:,:) = phytot .* (PDF .* xdel .* ydel);
    phy_disc(j,:,:) = PHYsspdisc3D(depth(j),SamplingDays(j),:,:);
end
%............................................................................
%----------------------------------------------------------------------------

figure(fignum)
%............................................................................
subplot(1,2,1)
counter = 0;
hplot001 = zeros(counter,1);
hplot002 = zeros(counter,1);
hplot003 = zeros(counter,1);
hplot004 = zeros(counter,1);
% for jday = SamplingDays(:)'
for jday = 1:3
    counter =  counter + 1;
    % phy_cont_dayj = phy_cont(jday,:);
    phy_cont_dayj = sum(phy_cont(jday,:,:),3); % Le Gland, 05/09/2019
    phy_disc_dayj = sum(phy_disc(jday,:,:),3); % Le Gland, 05/09/2019
    hplot001(counter) = plot(xrng(xIndex),phy_cont_dayj(xIndex),'k-');
    hold on
    hplot002(counter) = plot(xrng(xIndex),phy_cont_dayj(xIndex),'.');
    hold on
    hplot003(counter) = plot(xrng(xIndex),phy_disc_dayj(xIndex),'k--');
    hold on
    hplot004(counter) = plot(xrng(xIndex),phy_disc_dayj(xIndex),'.');
    color_dayj = cmap(counter,:);
    set(hplot001(counter),'Color',color_dayj)
    set(hplot002(counter),'Marker','p','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color_dayj,'MarkerSize',Markersize*2)
    set(hplot003(counter),'Color',color_dayj*0.6)
    set(hplot004(counter),'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color_dayj*0.6,'MarkerSize',Markersize)
end
hold off
grid on 
% legend(hplot002,num2str(SamplingDays(:)),'Location','NorthWest');
legend(hplot004,['z =  55m, day = 135';'z =   5m, day = 225';'z = 105m, day = 300'],'Location','NorthWest');
%%set(gca,'Xtick',xrng,'XtickLabel',exp(xrng))
set(gca,'Ylim',[ymin ymax])
set(gca,'Fontsize',Fontsize)
xlabel('log(Kd)')
ylabel('Biomass')
% title('Phytoplankton Kd distribution -- Cont','Fontsize',Fontsize + 2)
title('Phytoplankton Kd distribution','Fontsize',Fontsize + 2)
%............................................................................
% subplot(2,2,2)
% counter = 0;
% for jday = SamplingDays(:)'
%     counter =  counter + 1;
%     %phy_disc_dayj = phy_disc(jday,:);
%     phy_disc_dayj = sum(phy_disc(jday,:,:),3); % Le Gland, 05/09/2019
%     hplot001 = plot(xrng(xIndex),phy_disc_dayj(xIndex),'k-');
%     hold on
%     hplot002 = plot(xrng(xIndex),phy_disc_dayj(xIndex),'.');
%     hold on
%     color_dayj = cmap(counter,:);
%     set(hplot001,'Color',color_dayj)
%     set(hplot002,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color_dayj,'MarkerSize',Markersize) 
% end
% hold off
% grid on 
% set(gca,'Ylim',[ymin ymax])
% set(gca,'Fontsize',Fontsize)
% xlabel('log(size)')
% ylabel('Biomass')
% title('Phytoplankton size distribution -- Disc','Fontsize',Fontsize + 2)

%............................................................................
% Plot the distribution of Topt (Le Gland, 05/09/2019)

%subplot(2,2,3)
subplot(1,2,2)
counter = 0;
% for jday = SamplingDays(:)'
for jday = 1:3
    counter =  counter + 1;
    phy_cont_dayj = squeeze(sum(phy_cont(jday,:,:),2));
    phy_disc_dayj = squeeze(sum(phy_disc(jday,:,:),2));
    hplot001(counter) = plot(yrng(yIndex),phy_cont_dayj(yIndex),'k-');
    hold on
    hplot002(counter) = plot(yrng(yIndex),phy_cont_dayj(yIndex),'.');
    hold on
    hplot003(counter) = plot(yrng(yIndex),phy_disc_dayj(yIndex),'k--');
    hold on
    hplot004(counter) = plot(yrng(yIndex),phy_disc_dayj(yIndex),'.');
    color_dayj = cmap(counter,:);
    set(hplot001(counter),'Color',color_dayj)
    set(hplot002(counter),'Marker','p','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color_dayj,'MarkerSize',Markersize*2)
    set(hplot003(counter),'Color',color_dayj*0.6)
    set(hplot004(counter),'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color_dayj*0.6,'MarkerSize',Markersize)
end
hold off
grid on 
%legend(hplot002,num2str(SamplingDays(:)),'Location','NorthWest');
legend(hplot004,['z =  55m, day = 135';'z =   5m, day = 225';'z = 105m, day = 300'],'Location','NorthWest');
%%set(gca,'Xtick',xrng,'XtickLabel',exp(xrng))
set(gca,'Ylim',[ymin ymax])
set(gca,'Fontsize',Fontsize)
xlabel('Topt')
ylabel('Biomass')
% title('Phytoplankton optimal temperature distribution -- Cont','Fontsize',Fontsize + 2)
title('Phytoplankton optimal temperature distribution','Fontsize',Fontsize + 2)
%............................................................................
% subplot(2,2,4)
% counter = 0;
% for jday = SamplingDays(:)'
%     counter =  counter + 1;
%     phy_disc_dayj = squeeze(sum(phy_disc(jday,:,:),2));
%     hplot001 = plot(yrng(yIndex),phy_disc_dayj(yIndex),'k-');
%     hold on
%     hplot002 = plot(yrng(yIndex),phy_disc_dayj(yIndex),'.');
%     hold on
%     color_dayj = cmap(counter,:);
%     set(hplot001,'Color',color_dayj)
%     set(hplot002,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color_dayj,'MarkerSize',Markersize) 
% end
% hold off
% grid on 
% set(gca,'Ylim',[ymin ymax])
% set(gca,'Fontsize',Fontsize)
% xlabel('Topt')
% ylabel('Biomass')
% title('Phytoplankton optimal temperature distribution -- Disc','Fontsize',Fontsize + 2)
%............................................................................
%============================================================================
%****************************************************************************
return

    
