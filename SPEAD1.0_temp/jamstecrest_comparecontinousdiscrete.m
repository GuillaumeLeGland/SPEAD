
%===================================================================================
% $$$ %...................................................................................
% $$$ xtime = todedotday;
% $$$ %...................................................................................
% $$$ PHYTcontData = PHYTsspcont(:);
% $$$ PHYTdiscData = PHYTsspdisc(:); 
% $$$ ZOOcontData = ZOOsspcont(:);
% $$$ ZOOdiscData = ZOOsspdisc(:); 
% $$$ DINcontData = DINsspcont(:);
% $$$ DINdiscData = DINsspdisc(:); 
% $$$ PONcontData = PONsspcont(:);
% $$$ PONdiscData = PONsspdisc(:); 
% $$$ %...................................................................................
% $$$ FPHYTcontData = FPHYTsspcont(:);
% $$$ FPHYTdiscData = FPHYTsspdisc(:); 
% $$$ FZOOcontData = FZOOsspcont(:);
% $$$ FZOOdiscData = FZOOsspdisc(:); 
% $$$ FDINcontData = FDINsspcont(:);
% $$$ FDINdiscData = FDINsspdisc(:); 
% $$$ FPONcontData = FPONsspcont(:);
% $$$ FPONdiscData = FPONsspdisc(:); 
% $$$ %...................................................................................
% $$$ UXData = UXssp(:);
% $$$ GXData = GXssp(:);
% $$$ uphytotData = uphytotssp(:);
% $$$ gphytotData = gphytotssp(:); 
% $$$ %...................................................................................
% $$$ logESDphyAveContData = logESDphysspAveCont(:);  
% $$$ logESDphyAveDiscData = logESDphysspAveDisc(:);  
% $$$ logESDphyStdContData = logESDphysspStdCont(:); 
% $$$ logESDphyStdDiscData = logESDphysspStdDisc(:); 
% $$$ ESDphyAveContData = ESDphysspAveCont(:); 
% $$$ ESDphyAveDiscData = ESDphysspAveDisc(:);
% $$$ ESDphyCVcontData = ESDphysspCVcont(:); 
% $$$ ESDphyCVdiscData = ESDphysspCVdisc(:);
%...................................................................................
%===================================================================================
% $$$ %...................................................................................
% $$$ jdepth = 1; 
% $$$ %...................................................................................
% $$$ xtime = todedotday;
% $$$ %...................................................................................
% $$$ PHYTcontData = PHYTsspcont(jdepth,:);
% $$$ PHYTdiscData = PHYTsspdisc(jdepth,:); 
% $$$ ZOOcontData = ZOOsspcont(jdepth,:);
% $$$ ZOOdiscData = ZOOsspdisc(jdepth,:); 
% $$$ DINcontData = DINsspcont(jdepth,:);
% $$$ DINdiscData = DINsspdisc(jdepth,:); 
% $$$ PONcontData = PONsspcont(jdepth,:);
% $$$ PONdiscData = PONsspdisc(jdepth,:); 
% $$$ %...................................................................................
% $$$ FPHYTcontData = FPHYTsspcont(jdepth,:);
% $$$ FPHYTdiscData = FPHYTsspdisc(jdepth,:); 
% $$$ FZOOcontData = FZOOsspcont(jdepth,:);
% $$$ FZOOdiscData = FZOOsspdisc(jdepth,:); 
% $$$ FDINcontData = FDINsspcont(jdepth,:);
% $$$ FDINdiscData = FDINsspdisc(jdepth,:); 
% $$$ FPONcontData = FPONsspcont(jdepth,:);
% $$$ FPONdiscData = FPONsspdisc(jdepth,:); 
% $$$ %...................................................................................
% $$$ UXData = UXssp(jdepth,:);
% $$$ GXData = GXssp(jdepth,:);
% $$$ uphytotData = uphytotssp(jdepth,:);
% $$$ gphytotData = gphytotssp(jdepth,:); 
% $$$ %...................................................................................
% $$$ logESDphyAveContData = logESDphysspAveCont(jdepth,:);  
% $$$ logESDphyAveDiscData = logESDphysspAveDisc(jdepth,:);  
% $$$ logESDphyStdContData = logESDphysspStdCont(jdepth,:); 
% $$$ logESDphyStdDiscData = logESDphysspStdDisc(jdepth,:); 
% $$$ ESDphyAveContData = ESDphysspAveCont(jdepth,:); 
% $$$ ESDphyAveDiscData = ESDphysspAveDisc(jdepth,:);
% $$$ ESDphyCVcontData = ESDphysspCVcont(jdepth,:); 
% $$$ ESDphyCVdiscData = ESDphysspCVdisc(jdepth,:);
% $$$ %...................................................................................
%===================================================================================
%...................................................................................
%%jdepth = 1; 
%%jdepth = 5; 
%%jdepth = 8; 
jdepth = ndepths; 
%...................................................................................
if length(Jstep) == nsteps 
    xtime = todedotout; %For all years. 
elseif length(Jstep) == (nsteps/nyear)
    xtime = [deltat:deltat:ndays]; %For last year only. 
end 
%...................................................................................
PHYTcontData = PHYThdpcont(jdepth,:);
PHYTdiscData = PHYThdpdisc(jdepth,:); 
ZOOcontData = ZOOhdpcont(jdepth,:);
ZOOdiscData = ZOOhdpdisc(jdepth,:); 
DINcontData = DINhdpcont(jdepth,:);
DINdiscData = DINhdpdisc(jdepth,:); 
PONcontData = PONhdpcont(jdepth,:);
PONdiscData = PONhdpdisc(jdepth,:); 
%...................................................................................
FPHYTcontData = FPHYThdpcont(jdepth,:);
FPHYTdiscData = FPHYThdpdisc(jdepth,:); 
FZOOcontData = FZOOhdpcont(jdepth,:);
FZOOdiscData = FZOOhdpdisc(jdepth,:); 
FDINcontData = FDINhdpcont(jdepth,:);
FDINdiscData = FDINhdpdisc(jdepth,:); 
FPONcontData = FPONhdpcont(jdepth,:);
FPONdiscData = FPONhdpdisc(jdepth,:); 
%...................................................................................
UXData = UXhdp(jdepth,:);
GXData = GXhdp(jdepth,:);
%...................................................................................
uphytotData = uphytothdp(jdepth,:);
gphytotData = gphytothdp(jdepth,:); 
%...................................................................................
logESDphyAveContData = logESDphyhdpAveCont(jdepth,:);  
logESDphyAveDiscData = logESDphyhdpAveDisc(jdepth,:);  
logESDphyStdContData = logESDphyhdpStdCont(jdepth,:); 
logESDphyStdDiscData = logESDphyhdpStdDisc(jdepth,:); 
ESDphyAveContData = ESDphyhdpAveCont(jdepth,:); 
ESDphyAveDiscData = ESDphyhdpAveDisc(jdepth,:);
ESDphyCVcontData = ESDphyhdpCVcont(jdepth,:); 
ESDphyCVdiscData = ESDphyhdpCVdisc(jdepth,:);
%...................................................................................
%===================================================================================
% $$$ figure(300)
% $$$ %...................................................................................
% $$$ subplot(2,2,1) 
% $$$ plot(PHYTcontData,PHYTdiscData,'.')
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 2.0, 0 2.0])
% $$$ xlabel('PHY continous') 
% $$$ ylabel('PHY discrete') 
% $$$ title('PHY') 
% $$$ grid on
% $$$ subplot(2,2,2) 
% $$$ plot(ZOOcontData,ZOOdiscData,'.')
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 2.0, 0 2.0])
% $$$ xlabel('ZOO continous') 
% $$$ ylabel('ZOO discrete') 
% $$$ title('ZOO') 
% $$$ grid on
% $$$ subplot(2,2,3) 
% $$$ plot(DINcontData,DINdiscData,'.')
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 1.0, 0 1.0])
% $$$ xlabel('DIN continous') 
% $$$ ylabel('DIN discrete') 
% $$$ title('DIN') 
% $$$ grid on
% $$$ subplot(2,2,4) 
% $$$ plot(PONcontData,PONdiscData,'.')
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 5.0, 0 5.0])
% $$$ xlabel('PON continous') 
% $$$ ylabel('PON discrete') 
% $$$ title('PON') 
% $$$ grid on
% $$$ %...................................................................................
% $$$ %===================================================================================
% $$$ figure(310)
% $$$ %...................................................................................
% $$$ subplot(2,2,1) 
% $$$ %%plot(FPHYTcontBisData,FPHYTdiscData,'.')
% $$$ plot(FPHYTcontData,FPHYTdiscData,'.')
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 2.0, 0 2.0])
% $$$ xlabel('FPHY continous') 
% $$$ ylabel('FPHY discrete') 
% $$$ title('FPHY') 
% $$$ grid on
% $$$ subplot(2,2,2) 
% $$$ plot(FZOOcontData,FZOOdiscData,'.')
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 2.0, 0 2.0])
% $$$ xlabel('FZOO continous') 
% $$$ ylabel('FZOO discrete') 
% $$$ title('FZOO') 
% $$$ grid on
% $$$ subplot(2,2,3) 
% $$$ plot(FDINcontData,FDINdiscData,'.')
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 1.0, 0 1.0])
% $$$ xlabel('FDIN continous') 
% $$$ ylabel('FDIN discrete') 
% $$$ title('FDIN') 
% $$$ grid on
% $$$ subplot(2,2,4) 
% $$$ plot(FPONcontData,FPONdiscData,'.')
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 1.0, 0 1.0])
% $$$ xlabel('FPON continous') 
% $$$ ylabel('FPON discrete') 
% $$$ title('FPON') 
% $$$ grid on
% $$$ %...................................................................................
% $$$ %===================================================================================
% $$$ figure(320)
% $$$ %...................................................................................
% $$$ subplot(2,2,1) 
% $$$ %%plot(FPHYTcontData./PHYTcontData,FPHYTdiscData./PHYTdiscData,'.')
% $$$ plot(UXData,uphytotData,'.') 
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 1.2, 0 1.2])
% $$$ xlabel('FPHY/PHY continous') 
% $$$ ylabel('FPHY/PHY discrete') 
% $$$ title('FPHY/PHY') 
% $$$ grid on
% $$$ subplot(2,2,2) 
% $$$ %%plot(FZOOcontData./PHYTcontData,FZOOdiscData./PHYTdiscData,'.')
% $$$ plot(GXData,gphytotData,'.') 
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 1.2, 0 1.2])
% $$$ xlabel('FZOO/PHY continous') 
% $$$ ylabel('FZOO/PHY discrete') 
% $$$ title('FZOO/PHY') 
% $$$ grid on
% $$$ subplot(2,2,3) 
% $$$ plot(FDINcontData./DINcontData,FDINdiscData./DINdiscData,'.')
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis([0 1.0, 0 1.0])
% $$$ xlabel('FDIN/DIN continous') 
% $$$ ylabel('FDIN/DIN discrete') 
% $$$ title('FDIN/DIN') 
% $$$ grid on
% $$$ subplot(2,2,4) 
% $$$ plot(FPONcontData./PONcontData,FPONdiscData./PONcontData,'.')
% $$$ hold on
% $$$ plot([0:0.1:5.0],[0:0.1:5.0],'k-')
% $$$ hold off 
% $$$ %%axis([0 1.0, 0 1.0])
% $$$ xlabel('FPON/PON continous') 
% $$$ ylabel('FPON/PON discrete') 
% $$$ title('FPON/PON') 
% $$$ grid on
% $$$ %...................................................................................
% $$$ %===================================================================================
% $$$ %...................................................................................
% $$$ figure(330)
% $$$ subplot(2,2,1)
% $$$ plot(logESDphyAveContData,logESDphyAveDiscData,'*')
% $$$ hold on
% $$$ plot([0:0.1:10.0],[0:0.1:10.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 3, 0 3]) 
% $$$ xlabel('logESD ave - continous') 
% $$$ ylabel('logESD ave - discrete') 
% $$$ title('logESD ave')
% $$$ grid on
% $$$ subplot(2,2,2)
% $$$ plot(logESDphyStdContData,logESDphyStdDiscData,'*')
% $$$ hold on
% $$$ plot([0:0.1:10.0],[0:0.1:10.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 1, 0 1]) 
% $$$ xlabel('logESD std - continous') 
% $$$ ylabel('logESD std - discrete') 
% $$$ title('logESD std')
% $$$ grid on
% $$$ subplot(2,2,3)
% $$$ plot(ESDphyAveContData,ESDphyAveDiscData,'*')
% $$$ hold on
% $$$ plot([0:0.1:10.0],[0:0.1:10.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 10, 0 10]) 
% $$$ xlabel('ESD ave - continous') 
% $$$ ylabel('ESD ave - discrete') 
% $$$ title('ESD ave')
% $$$ grid on
% $$$ subplot(2,2,4)
% $$$ plot(ESDphyCVcontData,ESDphyCVdiscData,'*')
% $$$ hold on
% $$$ plot([0:0.1:10.0],[0:0.1:10.0],'k-')
% $$$ hold off 
% $$$ %%axis[0 3, 0 3]) 
% $$$ xlabel('ESD ave - continous') 
% $$$ ylabel('ESD ave - discrete') 
% $$$ title('ESD ave')
% $$$ grid on
% $$$ %...................................................................................
%===================================================================================
figure(400)
%...................................................................................
subplot(2,2,1) 
plot(xtime,PHYTcontData,'.r-',xtime,PHYTdiscData,'.b-')
axis([0 max(xtime), 0 1.2*PHYmax]) 
legend('continous','discrete')
title('PHY') 
grid on
subplot(2,2,2) 
plot(xtime,ZOOcontData,'.r-',xtime,ZOOdiscData,'.b-')
axis([0 max(xtime), 0 1.2*ZOOmax]) 
legend('continous','discrete')
title('ZOO') 
grid on
subplot(2,2,3) 
plot(xtime,DINcontData,'.r-',xtime,DINdiscData,'.b-')
axis([0 max(xtime), 0 1.2*DINmax]) 
legend('continous','discrete')
title('DIN') 
grid on
subplot(2,2,4) 
plot(xtime,PONcontData,'.r-',xtime,PONdiscData,'.b-')
axis([0 max(xtime), 0 1.2*PONmax]) 
legend('continous','discrete')
title('PON') 
grid on
%...................................................................................
print('-dpng','-r300','jamstecrest_continous-vs-discrete-0D_fig001.png')
%...................................................................................
%===================================================================================
figure(410)
%...................................................................................
subplot(2,2,1) 
plot(xtime,FPHYTcontData,'.r-',xtime,FPHYTdiscData,'.b-')
legend('continous','discrete')
title('FPHY') 
grid on
subplot(2,2,2) 
plot(xtime,FZOOcontData,'.r-',xtime,FZOOdiscData,'.b-')
legend('continous','discrete')
title('FZOO') 
grid on
subplot(2,2,3) 
plot(xtime,FDINcontData,'.r-',xtime,FDINdiscData,'.b-')
legend('continous','discrete')
title('FDIN') 
grid on
subplot(2,2,4) 
plot(xtime,FPONcontData,'.r-',xtime,FPONdiscData,'.b-')
legend('continous','discrete')
title('FPON') 
grid on
%...................................................................................
print('-dpng','-r300','jamstecrest_continous-vs-discrete-0D_fig002.png')
%...................................................................................
%===================================================================================
figure(420)
%...................................................................................
subplot(2,2,1) 
plot(xtime,UXData,'.r-',xtime,uphytotData,'.b-') 
legend('continous','discrete')
title('FPHY/PHY') 
grid on
subplot(2,2,2) 
plot(xtime,GXData,'.r-',xtime,gphytotData,'.b-') 
legend('continous','discrete')
title('FZOO/PHY') 
grid on
subplot(2,2,3) 
plot(xtime,FDINcontData./DINcontData,'.r-',xtime,FDINdiscData./DINdiscData,'.b-')
legend('continous','discrete')
title('FDIN/DIN') 
grid on
subplot(2,2,4) 
plot(xtime,FPONcontData./PONcontData,'.r-',xtime,FPONdiscData./PONcontData,'.b-')
legend('continous','discrete')
title('FPON/PON') 
grid on
%...................................................................................
print('-dpng','-r300','jamstecrest_continous-vs-discrete-0D_fig003.png')
%...................................................................................
%===================================================================================
%...................................................................................
figure(430)
subplot(2,2,1)
plot(xtime,logESDphyAveContData,'.r-',xtime,logESDphyAveDiscData,'-b.')
axis([0 max(xtime), 0 log(ESDmax/2)]) 
legend('continous','discrete')
title('logESD ave')
grid on
subplot(2,2,2)
plot(xtime,logESDphyStdContData,'.r-',xtime,logESDphyStdDiscData,'-b.')
axis([0 max(xtime), 0 log(ESDmax/2)]) 
legend('continous','discrete')
title('logESD std')
grid on
subplot(2,2,3)
plot(xtime,ESDphyAveContData,'.r-',xtime,ESDphyAveDiscData,'-b.')
axis([0 max(xtime), 0 ESDmax/2]) 
legend('continous','discrete')
title('ESD ave')
grid on
subplot(2,2,4)
plot(xtime,ESDphyCVcontData,'.r-',xtime,ESDphyCVdiscData,'-b.')
axis([0 max(xtime), 0 ESDmax/2]) 
legend('continous','discrete')
title('ESD ave')
grid on
%...................................................................................
print('-dpng','-r300','jamstecrest_continous-vs-discrete-0D_fig004.png')
%...................................................................................
%===================================================================================

