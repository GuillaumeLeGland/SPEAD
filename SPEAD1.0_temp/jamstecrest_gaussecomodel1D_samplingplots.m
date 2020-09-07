%function [] = jamstecrest_gaussecomodel1D_samplingplots(phy_cont,phy_disc,xrng,ndays,fignum)
% Function changed by Le Gland (05/09/2019) to plot 2 traits
function [] = jamstecrest_gaussecomodel1D_samplingplots(phy_cont,phy_disc,xrng,yrng,ndays,fignum) 

%============================================================================
%............................................................................
%%SamplingDays = [1,[5:5:ndays/2]];
%%SamplingDays = [1,[20:20:ndays/2]];
%%SamplingDays = [1,[15:15:ndays/2]];
SamplingDays = [1,[15:15:2*ndays/3]];
%%SamplingDays = [1,[30:30:ndays/2]];
nsamplingdays = length(SamplingDays);
%............................................................................
cmap = [jet(nsamplingdays)];
%%cmap = [jet(nsamplingdays/2);flipud(jet(nsamplingdays/2))]; %Symmetric in time.
%............................................................................
Fontsize = [8];
Markersize = [5];
%............................................................................
%%Index = [2:2:64];
%Index = [1:1:64];
%Index = [1:1:size(xrng,2)]; % Le Gland (06/06/2019)
% Index differs for each trait (Le Gland, 12/09/2019)
xIndex = [1:1:size(xrng,2)];
yIndex = [1:1:size(yrng,2)];
%............................................................................
ymin = 0;
%ymax = 0.2;
ymax = 0.4;
%%ymax = 0.5;
%............................................................................
%============================================================================
figure(fignum)
%............................................................................
subplot(2,2,1)
counter = 0;
for jday = SamplingDays(:)'
    counter =  counter + 1;
    % phy_cont_dayj = phy_cont(jday,:);
    phy_cont_dayj = sum(phy_cont(jday,:,:),3); % Le Gland, 05/09/2019
    hplot001(counter) = plot(xrng(xIndex),phy_cont_dayj(xIndex),'k-');
    hold on
    hplot002(counter) = plot(xrng(xIndex),phy_cont_dayj(xIndex),'.');
    hold on
    color_dayj = cmap(counter,:);
    set(hplot001(counter),'Color',[color_dayj])
    set(hplot002(counter),'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color_dayj,'MarkerSize',Markersize) 
end
hold off
grid on 
legend(hplot002,num2str(SamplingDays(:)),'Location','NorthWest');
%%set(gca,'Xtick',xrng,'XtickLabel',exp(xrng))
set(gca,'Ylim',[ymin ymax])
set(gca,'Fontsize',Fontsize)
xlabel('log(size)')
ylabel('Biomass')
title('Phytoplankton size distribution -- Cont','Fontsize',[Fontsize + 2])
%............................................................................
subplot(2,2,2)
counter = 0;
for jday = SamplingDays(:)'
    counter =  counter + 1;
    %phy_disc_dayj = phy_disc(jday,:);
    phy_disc_dayj = sum(phy_disc(jday,:,:),3); % Le Gland, 05/09/2019
    hplot001 = plot(xrng(xIndex),phy_disc_dayj(xIndex),'k-');
    hold on
    hplot002 = plot(xrng(xIndex),phy_disc_dayj(xIndex),'.');
    hold on
    color_dayj = cmap(counter,:);
    set(hplot001,'Color',[color_dayj])
    set(hplot002,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color_dayj,'MarkerSize',Markersize) 
end
hold off
grid on 
set(gca,'Ylim',[ymin ymax])
set(gca,'Fontsize',Fontsize)
xlabel('log(size)')
ylabel('Biomass')
title('Phytoplankton size distribution -- Disc','Fontsize',[Fontsize + 2])

%............................................................................
% Plot the distribution of Topt (Le Gland, 05/09/2019)

subplot(2,2,3)
counter = 0;
for jday = SamplingDays(:)'
    counter =  counter + 1;
    phy_cont_dayj = squeeze(sum(phy_cont(jday,:,:),2));
    hplot001(counter) = plot(yrng(yIndex),phy_cont_dayj(yIndex),'k-');
    hold on
    hplot002(counter) = plot(yrng(yIndex),phy_cont_dayj(yIndex),'.');
    hold on
    color_dayj = cmap(counter,:);
    set(hplot001(counter),'Color',[color_dayj])
    set(hplot002(counter),'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color_dayj,'MarkerSize',Markersize) 
end
hold off
grid on 
legend(hplot002,num2str(SamplingDays(:)),'Location','NorthWest');
%%set(gca,'Xtick',xrng,'XtickLabel',exp(xrng))
set(gca,'Ylim',[ymin ymax])
set(gca,'Fontsize',Fontsize)
xlabel('Topt')
ylabel('Biomass')
title('Phytoplankton optimal temperature distribution -- Cont','Fontsize',[Fontsize + 2])
%............................................................................
subplot(2,2,4)
counter = 0;
for jday = SamplingDays(:)'
    counter =  counter + 1;
    phy_disc_dayj = squeeze(sum(phy_disc(jday,:,:),2));
    hplot001 = plot(yrng(yIndex),phy_disc_dayj(yIndex),'k-');
    hold on
    hplot002 = plot(yrng(yIndex),phy_disc_dayj(yIndex),'.');
    hold on
    color_dayj = cmap(counter,:);
    set(hplot001,'Color',[color_dayj])
    set(hplot002,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',color_dayj,'MarkerSize',Markersize) 
end
hold off
grid on 
set(gca,'Ylim',[ymin ymax])
set(gca,'Fontsize',Fontsize)
xlabel('Topt')
ylabel('Biomass')
title('Phytoplankton optimal temperature distribution -- Disc','Fontsize',[Fontsize + 2])
%............................................................................
%============================================================================
%****************************************************************************
return

    
