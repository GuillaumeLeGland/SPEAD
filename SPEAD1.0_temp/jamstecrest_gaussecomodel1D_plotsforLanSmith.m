function [] = jamstecrest_gaussecomodel1D_plotsforLanSmith(XAVEsspcont,XSTDsspcont,PHYsspcont,ZOOsspcont,DINsspcont,PONsspcont,UXssp,GXssp)
global galfa gbeta
global t0 deltat ndays nyear tmax 
global zdepths ndepths deltaz

%===================================================================================
%...................................................................................
logESDave = XAVEsspcont;
logESDstd = XSTDsspcont;
%...................................................................................
PHY = PHYsspcont;
ZOO = ZOOsspcont;
DIN = DINsspcont;
PON = PONsspcont;
%...................................................................................
NTOT = PHY + ZOO + DIN + PON;
%...................................................................................
UX = UXssp;
GX = GXssp;
%...................................................................................
%%galfaTarget = 2.0;
galfaTarget = 1.6;
%...................................................................................
distgalfa = (galfa - galfaTarget).^2; 
jdepth = find(distgalfa == min(distgalfa))
%...................................................................................
j001 = jdepth; %KTW yes.
jend = ndepths; %KTW not.
%...................................................................................
%===================================================================================
%...................................................................................
figure(1001)
subplot(2,2,1)
hplot001 = plot([1:ndays],PHY(j001,:),'.r-',[1:ndays],PHY(jend,:),'.b-');
set(gca,'Xlim',[0 360],'Xtick',[0:90:360],'XtickLabel',[0:90:360])
set(gca,'Ylim',[0 3.0],'Ytick',[0:01:3.0],'YtickLabel',[0:01:3.0])
title('Phy')
grid on
subplot(2,2,2)
hplot002 = plot([1:ndays],ZOO(j001,:),'.r-',[1:ndays],ZOO(jend,:),'.b-');
set(gca,'Xlim',[0 360],'Xtick',[0:90:360],'XtickLabel',[0:90:360])
set(gca,'Ylim',[0 1.6],'Ytick',[0:.4:1.6],'YtickLabel',[0:.4:1.6])
legend(hplot002,'KTW yes','KTW not')
title('Zoo')
grid on
subplot(2,2,3)
hplot003 = plot([1:ndays],DIN(j001,:),'.r-',[1:ndays],DIN(jend,:),'.b-');
set(gca,'Xlim',[0 360],'Xtick',[0:90:360],'XtickLabel',[0:90:360])
set(gca,'Ylim',[0 0.4],'Ytick',[0:.1:0.4],'YtickLabel',[0:.1:0.4])
title('DIN')
grid on
subplot(2,2,4)
hplot004 = plot([1:ndays],PON(j001,:),'.r-',[1:ndays],PON(jend,:),'.b-');
set(gca,'Xlim',[0 360],'Xtick',[0:90:360],'XtickLabel',[0:90:360])
set(gca,'Ylim',[0 4.0],'Ytick',[0:01:4.0],'YtickLabel',[0:01:4.0])
title('PON')
grid on
%...................................................................................
print('-depsc','-r300','jamstecrest_gaussecomodel1D_ForLanSmith-fig001.eps')
print('-dpng ','-r300','jamstecrest_gaussecomodel1D_ForLanSmith-fig001.png')
%===================================================================================
%...................................................................................
figure(1002)
subplot(2,2,1)
hplot001 = plot([1:ndays],UX(j001,:),'.r-',[1:ndays],UX(jend,:),'.b-');
set(gca,'Xlim',[0 360],'Xtick',[0:90:360],'XtickLabel',[0:90:360])
set(gca,'Ylim',[0 1.0],'Ytick',[0:.2:1.0],'YtickLabel',[0:.2:1.0])
title('Phy growth -- Uptake rate [d-1]')
grid on
subplot(2,2,2)
hplot002 = plot([1:ndays],GX(j001,:),'.r-',[1:ndays],GX(jend,:),'.b-');
set(gca,'Xlim',[0 360],'Xtick',[0:90:360],'XtickLabel',[0:90:360])
set(gca,'Ylim',[0 1.0],'Ytick',[0:.2:1.0],'YtickLabel',[0:.2:1.0])
title('Phy mortality -- Grazing rate [d-1]')
legend(hplot002,'KTW yes','KTW not')
grid on
subplot(2,2,3)
hplot003 = plot([1:ndays],logESDave(j001,:),'.r-',[1:ndays],logESDave(jend,:),'.b-');
set(gca,'Xlim',[0 360],'Xtick',[0:90:360],'XtickLabel',[0:90:360])
set(gca,'Ylim',[0 0.8],'Ytick',[0:.2:0.8],'YtickLabel',[0:.2:0.8])
title('Phy size -- log(ESD) ave')
grid on
subplot(2,2,4)
hplot003 = plot([1:ndays],logESDstd(j001,:),'.r-',[1:ndays],logESDstd(jend,:),'.b-');
set(gca,'Xlim',[0 360],'Xtick',[0:90:360],'XtickLabel',[0:90:360])
set(gca,'Ylim',[0 1.0],'Ytick',[0:.2:1.0],'YtickLabel',[0:.2:1.0])
title('Phy diversity -- log(ESD) std')
grid on
%...................................................................................
print('-depsc','-r300','jamstecrest_gaussecomodel1D_ForLanSmith-fig002.eps')
print('-dpng ','-r300','jamstecrest_gaussecomodel1D_ForLanSmith-fig002.png')
%===================================================================================

return

