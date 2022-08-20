

%expected inputs: key2T, keyDisc, keyModelResol, keyKN, keyTOPT

function []=SPEAD_1D_plot()

%JUST CHECKING IF PLOTING WORKS OKAY:
A256 = peaks(256);
A128x256 = A256(1:2:256,:);
A = A128x256;
Amin = min(A(:));
Amax = max(A(:));
fignum = 1001;
[hcbar] = SPEAD_1D_subplotesting(A,Amin,Amax,fignum,mypackages);
pause(0.5)
close all

% Plot external forcings
fignum = 14;
SPEAD_1D_imagescforcings(itemp,iparz0,PAR2D,imld,iKZ,fignum,mypackages)

% Plot trade-offs
fignum = 15;
SPEAD_1D_tradeoff(mup0,amup,[0.1,0.5,2.0],temp0,Q10a,18:4:30,fignum);

%........................................................................
figure(10)
plot(xaxis,fxtraitphy,'-b')
hold on
plot(xaxis,fxtraitphy,'r*')
hold off
grid on
%........................................................................
pause(0.5) 
close all 
pause(1)

figure(2020)
%...................................................................................
subplot(2,2,1)
plot(xaxis,fxtraitphy,'k-',xaxis,fxtraitphy,'b.')
hold on
plot(xless,0,'r*',xplus,0,'r*')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('log (size)')
ylabel('f (x)')
grid on
%...................................................................................
subplot(2,2,3)
plot(Kn,fxtraitphy,'k-',Kn,fxtraitphy,'b.')
hold on
plot(Knless,0,'r*',Knplus,0,'r*')
hold off
set(gca,'Xlim',[Knmin Knmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('log (half-sat)')
ylabel('f (x)')
grid on
%...................................................................................
pause(1)
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT STATISTICS OF BOTH CONTINUOUS AND DISCRETE MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%CONTINOUS:
if strcmp(key2T,'yes')
    %...................................................................................
    fignum = 1010;
    [hfig1010] = SPEAD_1D_imagescuptakerates(MUPsspcont,MUZsspcont,NTOTsspcont,ndepths,ndays,MUPmin,MUPmax,...
        MUZmin,MUZmax,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages);
    %...................................................................................
    fignum = 1020;
    [hfig1020] = SPEAD_1D_imagescNPZD(NTOTsspcont,CHLsspcont,PPsspcont,ZOOsspcont,DINsspcont,PONsspcont,ndepths,ndays,PHYmin,PHYmax,...
        ZOOmin,ZOOmax,DINmin,DINmax,PONmin,PONmax,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages);
    %...................................................................................
    % Compare model and observations
    fignum = 1022;
    [hfig1022] = SPEAD_gaussecomodel1D_imagescmodvsobs(12*(106/16)*PP_obs,CHL_obs,NO3_obs,PON_obs,12*(106/16)*PPsspcont,...
        CHLsspcont,DINsspcont,PONsspcont,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages);
    %...................................................................................
    fignum = 1030;
    [hfig1030] = SPEAD_1D_imagescstatistics(logESDphysspAveCont,logESDphysspStdCont,TOPTphysspAveCont,TOPTphysspStdCont,physspCorCont,PHYTsspcont,...
        ndepths,ndays,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,logESDaveMax,logESDaveMin,logESDstdMax,logESDstdMin,...
        TOPTaveMax,TOPTaveMin,TOPTstdMax,TOPTstdMin,CorrelationAbsMax,fignum,mypackages);
    %...................................................................................
end
%===================================================================================
%DISCRETE:
if strcmp(keyDisc,'yes')
    %...................................................................................
    fignum = 2010;
    [hfig2010] = SPEAD_1D_imagescuptakerates(MUPsspdisc,MUZsspdisc,NTOTsspdisc,ndepths,ndays,MUPmin,MUPmax,...
        MUZmin,MUZmax,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages);
    %...................................................................................
    fignum = 2020;
    [hfig2020] = SPEAD_1D_imagescNPZD(NTOTsspdisc,CHLsspdisc,PPsspdisc,ZOOsspdisc,DINsspdisc,PONsspdisc,ndepths,ndays,PHYmin,PHYmax,...
        ZOOmin,ZOOmax,DINmin,DINmax,PONmin,PONmax,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,fignum,mypackages);
    %...................................................................................
    fignum = 2030;
    [hfig2030] = SPEAD_1D_imagescstatistics(logESDphysspAveDisc,logESDphysspStdDisc,TOPTphysspAveDisc,TOPTphysspStdDisc,physspCorDisc,PHYTsspdisc,...
        ndepths,ndays,myXtickMarks,myXtickLabel,myYtickMarks,myYtickLabel,myYaxisLabel,logESDaveMax,logESDaveMin,logESDstdMax,logESDstdMin,...
        TOPTaveMax,TOPTaveMin,TOPTstdMax,TOPTstdMin,CorrelationAbsMax,fignum,mypackages);
    %...................................................................................
end
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTS OF TRAIT DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
if strcmp(key2T,'yes') && strcmp(keyDisc,'yes') && strcmp(keyModelResol,'1D')
    fignum = [70];
    SPEAD_1D_distribution(PHYTsspcont,logESDphysspAveCont,logESDphysspStdCont,TOPTphysspAveCont,TOPTphysspStdCont,physspCorCont,PHYsspdisc3D,xaxis,yaxis,itemp,DINsspdisc,fignum);
end
%...................................................................................
if strcmp(key2T,'yes') && strcmp(keyKN,'yes') && strcmp(keyTOPT,'yes')
    fignum = [80];
    SPEAD_gaussecomodel1D_surftraitplot(DINsspcont,DINsspcont_K,temp(:,1:360),...
    PHYTsspcont,PHYTsspcont_K,PHYTsspcont_T,XAVE_sspcont,XAVE_sspcont_K,XSTD_sspcont,XSTD_sspcont_K,...
    YAVE_sspcont,YAVE_sspcont_T,YSTD_sspcont,YSTD_sspcont_T,XYCOR_sspcont,fignum);
end
%...................................................................................
%===================================================================================
%...................................................................................
if strcmp(key2T,'yes') && strcmp(keyDisc,'yes')
    fignum = [24];
    SPEAD_1D_contvsdiscplot(PHYTsspdisc,PHYTsspcont,logESDphysspAveDisc,logESDphysspAveCont,...
    TOPTphysspAveDisc,TOPTphysspAveCont,logESDphysspStdDisc,logESDphysspStdCont,...
    TOPTphysspStdDisc,TOPTphysspStdCont,physspCorDisc,physspCorCont,ndepths,fignum);
end
%...................................................................................

end

