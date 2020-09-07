% Bell curves for cont and disc models at specific points (Le Gland,
% 18/11/2019)
function [] = jamstecrest_gaussecomodel1D_bellcurve(PHYTsspcont,logESDphysspAveCont,logESDphysspStdCont,TOPTphysspAveCont,TOPTphysspStdCont,physspCorCont,PHYsspdisc3D,xrng,yrng,fignum) 

%============================================================================
%............................................................................
Fontsize = 20;
Markersize = 5;
%............................................................................
xIndex = 1:1:length(xrng);
yIndex = 1:1:length(yrng);
%............................................................................
ymin = 0;
ymax = 0.16;
%............................................................................

%============================================================================

%============================================================================
%............................................................................
%Choose locations here (Le Gland, 19/11/2019)
depth = [6,1,11];
SamplingDays = [135,225,300];
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

    
