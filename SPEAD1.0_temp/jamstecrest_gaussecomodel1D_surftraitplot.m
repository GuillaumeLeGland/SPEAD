% Comparison between 1-trait and 2-traits models (Le Gland, 28/11/2019)
function [] = jamstecrest_gaussecomodel1D_surftraitplot(DINsspcont,DINsspcont_K,temp,...
PHYTsspcont,PHYTsspcont_K,PHYTsspcont_T,XAVE_sspcont,XAVE_sspcont_K,XSTD_sspcont,XSTD_sspcont_K,...
YAVE_sspcont,YAVE_sspcont_T,YSTD_sspcont,YSTD_sspcont_T,XYCOR_sspcont,myXtickMarks,myXtickLabel,fignum) 

%============================================================================
%............................................................................
% Fontsize = 20;
% Markersize = 5;
%............................................................................
% ymin = 0;
% ymax = 0.16;
jk = 1; % Choose depth level here 
%............................................................................

%============================================================================

%============================================================================
%............................................................................
figure(fignum)
%............................................................................
hplot = subplot(2,3,1);
hold off
plot(PHYTsspcont(jk,:),'g.')
hold on
plot(PHYTsspcont_K(jk,:),'b.')
plot(PHYTsspcont_T(jk,:),'r.')
set(hplot,'Xlim',[1 360  ],'Xtick',[0:30:360],'XtickLabel',num2str([0:30:360]'))
set(gca,'Ylim',[0.0 0.25])
xlabel('time [days]')
ylabel('phy conc [mmol m^{-3}]')
%legend({['2-trait model ';'Half-sat model';'T_{opt} model ']},'FontSize',12,'Location','NorthWest','Orientation','Horizontal');
grid on
%............................................................................
hplot = subplot(2,3,4);
hold off
plot(XYCOR_sspcont(jk,:),'g.')
hold on
plot([1,360],[0,0],'k-','linewidth',3)
set(hplot,'Xlim',[1 360  ],'Xtick',[0:30:360],'XtickLabel',num2str([0:30:360]'))
set(gca,'Ylim',[-1 1],'Ytick',-1:0.25:1)
xlabel('time [days]')
ylabel('Trait correlation')
grid on
%............................................................................
hplot = subplot(2,3,2);
hold off
plot(XAVE_sspcont(jk,:),'g.')
hold on
plot(XAVE_sspcont_K(jk,:),'b.')
plot(log(DINsspcont(jk,:)),'color',[0 0.5 0],'linestyle','--','linewidth',3)
plot(log(DINsspcont_K(jk,:)),'color',[0 0 0.5],'linestyle','--','linewidth',3)
set(hplot,'Xlim',[1 360  ],'Xtick',[0:30:360],'XtickLabel',num2str([0:30:360]'))
set(gca,'Ylim',[-2 0.5])
xlabel('time [days]')
ylabel('logKn ave')
grid on
%............................................................................
hplot = subplot(2,3,5);
hold off
plot(XSTD_sspcont(jk,:),'g.')
hold on
plot(XSTD_sspcont_K(jk,:),'b.')
set(hplot,'Xlim',[1 360  ],'Xtick',[0:30:360],'XtickLabel',num2str([0:30:360]'))
%set(gca,'Ylim',[0 0.8])
set(gca,'Ylim',[0 1.1])
xlabel('time [days]')
ylabel('logKn std')
grid on
%............................................................................
hplot = subplot(2,3,3);
hold off
plot(YAVE_sspcont(jk,:),'g.')
hold on
plot(YAVE_sspcont_T(jk,:),'r.')
plot(temp(jk,:),'k--','linewidth',3)
set(hplot,'Xlim',[1 360  ],'Xtick',[0:30:360],'XtickLabel',num2str([0:30:360]'))
set(gca,'Ylim',[19 28])
xlabel('time [days]')
ylabel('Topt ave [degree]')
grid on
%............................................................................
hplot = subplot(2,3,6);
hold off
plot(YSTD_sspcont(jk,:),'g.')
hold on
plot(YSTD_sspcont_T(jk,:),'r.')
set(hplot,'Xlim',[1 360  ],'Xtick',[0:30:360],'XtickLabel',num2str([0:30:360]'))
set(gca,'Ylim',[0 2],'Ytick',0:0.25:2)
xlabel('time [days]')
ylabel('Topt std [degree]')
grid on
%............................................................................
return

    
