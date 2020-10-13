function [] = SPEAD_1D_contvsdiscplot(PHYTsspdisc,PHYTsspcont,...
logESDphysspAveDisc,logESDphysspAveCont,TOPTphysspAveDisc,TOPTphysspAveCont,...
logESDphysspStdDisc,logESDphysspStdCont,TOPTphysspStdDisc,TOPTphysspStdCont,...
physspCorDisc,physspCorCont,ndepths,fignum)

% $$$ J = [1]; %surface only.
J = [1:ndepths]; %all depths.
figure(fignum)
%...................................................................................
% Phytoplankton disc/cont scatter plot
subplot (2,3,1)
plot(PHYTsspdisc(1:2:15,3:5:88),PHYTsspcont(1:2:15,3:5:88),'b*')
set(gca,'Xlim',[0 +0.3])
set(gca,'Ylim',[0 +0.3])
hold on
plot(PHYTsspdisc(1:2:15,93:5:178),PHYTsspcont(1:2:15,93:5:178),'g*')
plot(PHYTsspdisc(1:2:15,183:5:268),PHYTsspcont(1:2:15,183:5:268),'r*')
plot(PHYTsspdisc(1:2:15,273:5:358),PHYTsspcont(1:2:15,273:5:358),'k*')
plot([0,2],[0,2],'k-')
xlabel('phy conc -- Discrete')
ylabel('phy conc -- Continuous')
title('Phytoplankton concentration -- All depths')
R_phyt = corr(PHYTsspdisc(:),PHYTsspcont(:));
R2_phyt = R_phyt.^2;
annotation('textbox', [0.135, 0.82, 0.1, 0.1], 'String', ['R^2 = ', num2str(0.001*round(1000*R2_phyt))], 'FontSize', 15);
grid on
%...................................................................................
subplot(2,3,2)
plot(logESDphysspAveDisc(1:2:15,3:5:88),logESDphysspAveCont(1:2:15,3:5:88),'b*')
set(gca,'Xlim',[-2.5 +1.5])
set(gca,'Ylim',[-2.5 +1.5])
hold on
plot(logESDphysspAveDisc(1:2:15,93:5:178),logESDphysspAveCont(1:2:15,93:5:178),'g*')
plot(logESDphysspAveDisc(1:2:15,183:5:268),logESDphysspAveCont(1:2:15,183:5:268),'r*')
plot(logESDphysspAveDisc(1:2:15,273:5:358),logESDphysspAveCont(1:2:15,273:5:358),'k*')
plot([-2.5,1.5],[-2.5,1.5],'k-')
xlabel('ave logsize -- Discrete')
ylabel('ave logsize -- Continuous')
title('Mean size -- All depths')
R_logESDave = corr(logESDphysspAveDisc(:),logESDphysspAveCont(:));
R2_logESDave = R_logESDave.^2;
annotation('textbox', [0.415, 0.82, 0.1, 0.1], 'String', ['R^2 = ', num2str(0.001*round(1000*R2_logESDave))], 'FontSize', 15);
grid on
%...................................................................................
subplot(2,3,5)
plot(logESDphysspStdDisc(1:2:15,3:5:88),logESDphysspStdCont(1:2:15,3:5:88),'b*')
set(gca,'Xlim',[0.0 1.2])
set(gca,'Ylim',[0.0 1.2])
hold on
plot(logESDphysspStdDisc(1:2:15,93:5:178),logESDphysspStdCont(1:2:15,93:5:178),'g*')
plot(logESDphysspStdDisc(1:2:15,183:5:268),logESDphysspStdCont(1:2:15,183:5:268),'r*')
plot(logESDphysspStdDisc(1:2:15,273:5:358),logESDphysspStdCont(1:2:15,273:5:358),'k*')
plot([0.0,1.2],[0.0,1.2],'k-')
xlabel('std logsize -- Discrete')
ylabel('std logsize -- Continuous')
title('Standard Deviation logsize -- All depths')
R_logESDstd = corr(logESDphysspStdDisc(:),logESDphysspStdCont(:));
R2_logESDstd = R_logESDstd.^2;
annotation('textbox', [0.415, 0.345, 0.1, 0.1], 'String', ['R^2 = ', num2str(0.001*round(1000*R2_logESDstd))], 'FontSize', 15);
grid on

subplot(2,3,4)
plot(physspCorDisc(1:2:15,3:5:88),physspCorCont(1:2:15,3:5:88),'b*')
set(gca,'Xlim',[-1 1],'Xtick',-1:0.25:1)
set(gca,'Ylim',[-1 1],'Ytick',-1:0.25:1)
hold on
plot(physspCorDisc(1:2:15,93:5:178),physspCorCont(1:2:15,93:5:178),'g*')
plot(physspCorDisc(1:2:15,183:5:268),physspCorCont(1:2:15,183:5:268),'r*')
plot(physspCorDisc(1:2:15,273:5:358),physspCorCont(1:2:15,273:5:358),'k*')
plot([-1,1],[-1,1],'k-')
xlabel('Correlation -- Discrete')
ylabel('Correlation -- Continuous')
title('Correlation -- All depths')
R_Cor = corr(physspCorDisc(:),physspCorCont(:));
R2_Cor = R_Cor.^2;
annotation('textbox', [0.135, 0.345, 0.1, 0.1], 'String', ['R^2 = ', num2str(0.001*round(1000*R2_Cor))], 'FontSize', 15);
grid on
%...................................................................................
subplot(2,3,3)
plot(TOPTphysspAveDisc(1:2:15,3:5:88),TOPTphysspAveCont(1:2:15,3:5:88),'b*')
set(gca,'Xlim',[18 30])
set(gca,'Ylim',[18 30])
hold on
plot(TOPTphysspAveDisc(1:2:15,93:5:178),TOPTphysspAveCont(1:2:15,93:5:178),'g*')
plot(TOPTphysspAveDisc(1:2:15,183:5:268),TOPTphysspAveCont(1:2:15,183:5:268),'r*')
plot(TOPTphysspAveDisc(1:2:15,273:5:358),TOPTphysspAveCont(1:2:15,273:5:358),'k*')
plot([18,30],[18,30],'k-')
xlabel('ave Topt -- Discrete')
ylabel('ave Topt -- Continuous')
title('Mean Topt -- All depths')
R_TOPTave = corr(TOPTphysspAveDisc(:),TOPTphysspAveCont(:));
R2_TOPTave = R_TOPTave.^2;
annotation('textbox', [0.695, 0.82, 0.1, 0.1], 'String', ['R^2 = ', num2str(0.001*round(1000*R2_TOPTave))], 'FontSize', 15);
grid on
%...................................................................................
subplot(2,3,6)
plot(TOPTphysspStdDisc(1:2:15,3:5:88),TOPTphysspStdCont(1:2:15,3:5:88),'b*')
set(gca,'Xlim',[0.0 4.0])
set(gca,'Ylim',[0.0 4.0])
hold on
plot(TOPTphysspStdDisc(1:2:15,93:5:178),TOPTphysspStdCont(1:2:15,93:5:178),'g*')
plot(TOPTphysspStdDisc(1:2:15,183:5:268),TOPTphysspStdCont(1:2:15,183:5:268),'r*')
plot(TOPTphysspStdDisc(1:2:15,273:5:358),TOPTphysspStdCont(1:2:15,273:5:358),'k*')
plot([0.0,4.0],[0.0,4.0],'k-')
xlabel('std Topt -- Discrete')
ylabel('std Topt -- Continuous')
title('Standard Deviation Topt -- All depths')
R_TOPTstd = corr(TOPTphysspStdDisc(:),TOPTphysspStdCont(:));
R2_TOPTstd = R_TOPTstd.^2;
annotation('textbox', [0.695, 0.345, 0.1, 0.1], 'String', ['R^2 = ', num2str(0.001*round(1000*R2_TOPTstd))], 'FontSize', 15);
grid on