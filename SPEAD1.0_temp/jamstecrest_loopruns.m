%********************************************************************
%====================================================================
%....................................................................
close all
clear all
format short g 
%....................................................................
%====================================================================
%....................................................................
OnePath001 = '~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYFUNCTIONS/';
GenPath001 = genpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYTOOLBOX/');
%....................................................................
%%GenPath002 = genpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/TOOLBOX-MEX/NETCDF-MATLAB6p5/');
GenPath002 = genpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/TOOLBOX-MEX/NETCDF-MATLAB2011a-64bit/');
%....................................................................
addpath(OnePath001)
addpath(GenPath001)
addpath(GenPath002)
%....................................................................
%====================================================================
%********************************************************************
%....................................................................
global myYtickMarks myYtickLabel myYaxisLabel 
global myXtickMarks myXtickLabel myXaxisLabel 
%....................................................................
%====================================================================
%....................................................................
for jloop = 1:8 
    %================================================================
    %................................................................
    [SDINloop,...
     MUPcontloop,MUZcontloop,NTOTcontloop,...
     MUPdiscloop,MUZdiscloop,NTOTdiscloop,...
     PHYTcontloop,ZOOcontloop,DINcontloop,PONcontloop,...
     PHYTdiscloop,ZOOdiscloop,DINdiscloop,PONdiscloop,...
     FPHYTcontloop,FZOOcontloop,FDINcontloop,FPONcontloop,...
     FPHYTdiscloop,FZOOdiscloop,FDINdiscloop,FPONdiscloop,...
     logESDphyAveContloop,logESDphyStdContloop,ESDphyAveContloop,ESDphyCVcontloop,...
     logESDphyAveDiscloop,logESDphyStdDiscloop,ESDphyAveDiscloop,ESDphyCVdiscloop,...
     ShannonEntropyTheoreticalContloop,SindexContloop,EindexContloop,...
     ShannonEntropyTheoreticalDiscloop,SindexDiscloop,EindexDiscloop] = jamstecrest_gaussecomodel1D(jloop);
    %................................................................
    %================================================================
    %................................................................
    SDIN(:,:,jloop) = SDINloop;
    %................................................................
    %================================================================
    %................................................................
    MUPcont(:,:,jloop) = MUPcontloop;
    MUZcont(:,:,jloop) = MUZcontloop;
    NTOTcont(:,:,jloop) = NTOTcontloop;
    %................................................................
    PHYTcont(:,:,jloop) = PHYTcontloop;
    ZOOcont(:,:,jloop)  = ZOOcontloop;
    DINcont(:,:,jloop)  = DINcontloop;
    PONcont(:,:,jloop)  = PONcontloop;
    %................................................................
    FPHYTcont(:,:,jloop) = FPHYTcontloop;
    FZOOcont(:,:,jloop)  = FZOOcontloop;
    FDINcont(:,:,jloop)  = FDINcontloop;
    FPONcont(:,:,jloop)  = FPONcontloop;
    %................................................................
    logESDphyAveCont(:,:,jloop) = logESDphyAveContloop;
    logESDphyStdCont(:,:,jloop) = logESDphyStdContloop;
    ESDphyAveCont(:,:,jloop)    = ESDphyAveContloop;
    ESDphyCVcont(:,:,jloop)     = ESDphyCVcontloop;
    %................................................................
    ShannonEntropyTheoreticalCont(:,:,jloop) = ShannonEntropyTheoreticalContloop;
    SindexCont(:,:,jloop) = SindexContloop; 
    EindexCont(:,:,jloop) = EindexContloop; 
    %................................................................
    %================================================================
    %................................................................
    MUPdisc(:,:,jloop) = MUPdiscloop;
    MUZdisc(:,:,jloop) = MUZdiscloop;
    NTOTdisc(:,:,jloop) = NTOTdiscloop;
    %................................................................
    PHYTdisc(:,:,jloop) = PHYTdiscloop;
    ZOOdisc(:,:,jloop)  = ZOOdiscloop;
    DINdisc(:,:,jloop)  = DINdiscloop;
    PONdisc(:,:,jloop)  = PONdiscloop;
    %................................................................
    FPHYTdisc(:,:,jloop) = FPHYTdiscloop;
    FZOOdisc(:,:,jloop)  = FZOOdiscloop;
    FDINdisc(:,:,jloop)  = FDINdiscloop;
    FPONdisc(:,:,jloop)  = FPONdiscloop;
    %................................................................
    logESDphyAveDisc(:,:,jloop) = logESDphyAveDiscloop;
    logESDphyStdDisc(:,:,jloop) = logESDphyStdDiscloop;
    ESDphyAveDisc(:,:,jloop)    = ESDphyAveDiscloop;
    ESDphyCVdisc(:,:,jloop)     = ESDphyCVdiscloop;
    %................................................................
    ShannonEntropyTheoreticalDisc(:,:,jloop) = ShannonEntropyTheoreticalDiscloop;
    SindexDisc(:,:,jloop) = SindexDiscloop; 
    EindexDisc(:,:,jloop) = EindexDiscloop; 
    %................................................................
    %================================================================
end
%====================================================================
%??????????????????????????????????????
% $$$ MUZcont = MUPcont;
% $$$ MUZdisc = MUPdisc;
%??????????????????????????????????????
%....................................................................
clear *Path* 
clear *loop 
%....................................................................
%====================================================================
%....................................................................
% $$$ save('jamstecrest-loopruns-v000_DINflux-press_KTW.mat') 
%....................................................................
save('jamstecrest-loopruns-v000_DINflux-pulse_KTW.mat') 
%....................................................................
%====================================================================
*********************************************************************
return
%====================================================================
%....................................................................
load /home/vallina/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MODELLING/JAMSTEC/CODE-WORK/CODE-v055/MAT/galfa.txt
%....................................................................
% $$$ load /home/vallina/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MODELLING/JAMSTEC/CODE-WORK/CODE-v055/MAT/jamstecrest-loopruns-v055_DINflux-press_KTW.mat
%....................................................................
load /home/vallina/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MODELLING/JAMSTEC/CODE-WORK/CODE-v055/MAT/jamstecrest-loopruns-v055_DINflux-pulse_KTW.mat
%....................................................................
%====================================================================
%....................................................................
EindexChapCont = exp(ShannonEntropyTheoreticalCont);
%....................................................................
EindexChapDisc = exp(ShannonEntropyTheoreticalDisc);
%....................................................................
%====================================================================
%....................................................................
[aveSDIN,aveNTOTcont,aveMUPcont,aveMUZcont,aveFPHYTcont,aveFZOOcont,avePHYTcont,aveZOOcont,aveDINcont,avePONcont,avelogESDphyAveCont,avelogESDphyStdCont,aveESDphyAveCont,aveESDphyCVcont,aveSindexCont,aveEindexCont,aveEindexChapCont] = jamstecrest_loopruns_longtermaverage(SDIN,NTOTcont,MUPcont,MUZcont,FPHYTcont,FZOOcont,PHYTcont,ZOOcont,DINcont,PONcont,logESDphyAveCont,logESDphyStdCont,ESDphyAveCont,ESDphyCVcont,SindexCont,EindexCont,EindexChapCont);
%....................................................................
[aveSDIN,aveNTOTdisc,aveMUPdisc,aveMUZdisc,aveFPHYTdisc,aveFZOOdisc,avePHYTdisc,aveZOOdisc,aveDINdisc,avePONdisc,avelogESDphyAveDisc,avelogESDphyStdDisc,aveESDphyAveDisc,aveESDphyCVdisc,aveSindexDisc,aveEindexDisc,aveEindexChapDisc] = jamstecrest_loopruns_longtermaverage(SDIN,NTOTdisc,MUPdisc,MUZdisc,FPHYTdisc,FZOOdisc,PHYTdisc,ZOOdisc,DINdisc,PONdisc,logESDphyAveDisc,logESDphyStdDisc,ESDphyAveDisc,ESDphyCVdisc,SindexDisc,EindexDisc,EindexChapDisc);
%....................................................................
%====================================================================
%....................................................................
[aveSDINI] = jamstecrest_loopruns_interpolate(aveSDIN);
%....................................................................
[aveMUPcontI] = jamstecrest_loopruns_interpolate(aveMUPcont);
[aveMUPdiscI] = jamstecrest_loopruns_interpolate(aveMUPdisc);
%....................................................................
[aveMUZcontI] = jamstecrest_loopruns_interpolate(aveMUZcont);
[aveMUZdiscI] = jamstecrest_loopruns_interpolate(aveMUZdisc);
%....................................................................
[aveNTOTcontI] = jamstecrest_loopruns_interpolate(aveNTOTcont);
[aveNTOTdiscI] = jamstecrest_loopruns_interpolate(aveNTOTdisc);
%....................................................................
[avePHYTcontI] = jamstecrest_loopruns_interpolate(avePHYTcont);
[avePHYTdiscI] = jamstecrest_loopruns_interpolate(avePHYTdisc);
%....................................................................
[aveZOOcontI] = jamstecrest_loopruns_interpolate(aveZOOcont);
[aveZOOdiscI] = jamstecrest_loopruns_interpolate(aveZOOdisc);
%....................................................................
[aveFPHYTcontI] = jamstecrest_loopruns_interpolate(aveFPHYTcont);
[aveFPHYTdiscI] = jamstecrest_loopruns_interpolate(aveFPHYTdisc);
%....................................................................
[aveFZOOcontI] = jamstecrest_loopruns_interpolate(aveFZOOcont);
[aveFZOOdiscI] = jamstecrest_loopruns_interpolate(aveFZOOdisc);
%....................................................................
[aveDINcontI] = jamstecrest_loopruns_interpolate(aveDINcont);
[aveDINdiscI] = jamstecrest_loopruns_interpolate(aveDINdisc);
%....................................................................
[avePONcontI] = jamstecrest_loopruns_interpolate(avePONcont);
[avePONdiscI] = jamstecrest_loopruns_interpolate(avePONdisc);
%....................................................................
[aveESDphyAveContI] = jamstecrest_loopruns_interpolate(aveESDphyAveCont);
[aveESDphyAveDiscI] = jamstecrest_loopruns_interpolate(aveESDphyAveDisc);
%....................................................................
[aveESDphyCVcontI] = jamstecrest_loopruns_interpolate(aveESDphyCVcont);
[aveESDphyCVdiscI] = jamstecrest_loopruns_interpolate(aveESDphyCVdisc);
%....................................................................
[avelogESDphyAveContI] = jamstecrest_loopruns_interpolate(avelogESDphyAveCont);
[avelogESDphyAveDiscI] = jamstecrest_loopruns_interpolate(avelogESDphyAveDisc);
%....................................................................
[avelogESDphyStdContI] = jamstecrest_loopruns_interpolate(avelogESDphyStdCont);
[avelogESDphyStdDiscI] = jamstecrest_loopruns_interpolate(avelogESDphyStdDisc);
%....................................................................
[aveSindexContI] = jamstecrest_loopruns_interpolate(aveSindexCont);
[aveSindexDiscI] = jamstecrest_loopruns_interpolate(aveSindexDisc);
%....................................................................
[aveEindexContI] = jamstecrest_loopruns_interpolate(aveEindexCont);
[aveEindexDiscI] = jamstecrest_loopruns_interpolate(aveEindexDisc);
%....................................................................
[aveEindexChapContI] = jamstecrest_loopruns_interpolate(aveEindexChapCont);
[aveEindexChapDiscI] = jamstecrest_loopruns_interpolate(aveEindexChapDisc);
%....................................................................
[msize,nsize] = size(aveSDINI);
%....................................................................
%====================================================================
%....................................................................
% $$$ delgalfa = (max(galfa) - min(galfa)) / (msize - 1);
% $$$ %....................................................................
% $$$ galfaIbis = fliplr([min(galfa):delgalfa:max(galfa)]); %WRONG!!!
% $$$ %....................................................................
%====================================================================
%....................................................................
z = [1:length(galfa)];
dz = (max(z) - min(z)) / (msize - 1);
zi = [min(z):dz:max(z)]';
%....................................................................
galfaI = interp1(z,galfa,zi);
%....................................................................
%====================================================================
%....................................................................
nbins = 4;
ndepths = msize;
dn = ndepths/nbins;
%....................................................................
J = [1,[dn:dn:ndepths]];
%....................................................................
galfaJ = galfaI(J);
galfaJround = ceil(galfaJ*10)/10;
%....................................................................
myYtickMarks = [1,[dn:dn:ndepths]];
myYtickLabel = num2str(galfaJround); %(for vertical gradient in KTW alfa parameter with 0D model)
myYaxisLabel = 'KTW alfa';
%....................................................................
%====================================================================
%....................................................................
DoublingSeries = [1,2,4,8,16,32,64,128,256,512];
%....................................................................
PulseFrequency = DoublingSeries(:);
%....................................................................
nfreqs = length(PulseFrequency);
%....................................................................
x = [1:nfreqs];
dx = (max(x) - min(x)) / (msize - 1);
xi = [min(x):dx:max(x)]';
%....................................................................
PulseFrequencyI = interp1(x,PulseFrequency,xi);
%....................................................................
%====================================================================
%....................................................................
nbins = 4;
ntimes = msize;
dn = ntimes/nbins;
%....................................................................
%%Jticks = [1,[dn:dn:ntimes]];
Jticks = [1,14,27,40];
%....................................................................
PulseFreqJ = PulseFrequencyI(Jticks);
PulseFreqJround = ceil(PulseFreqJ*10)/10;
%....................................................................
myXtickMarks = Jticks;
myXtickLabel = num2str(PulseFreqJround); %(for vertical gradient in KTW alfa parameter with 0D model)
myXaxisLabel = 'Pulse frequency (per year)';
%....................................................................
%====================================================================
%....................................................................
aveDiversityCont  = aveSindexCont;
aveDiversityDisc  = aveSindexDisc;
%....................................................................
% $$$ aveDiversityCont  = aveEindexCont;
% $$$ aveDiversityDisc  = aveEindexDisc;
%....................................................................
% $$$ aveDiversityCont  = aveEindexChapCont;
% $$$ aveDiversityDisc  = aveEindexChapDisc;
%....................................................................
%====================================================================
%....................................................................
aveDiversityContI = aveSindexContI;
aveDiversityDiscI = aveSindexDiscI;
%....................................................................
% $$$ aveDiversityContI = aveEindexContI;
% $$$ aveDiversityDiscI = aveEindexDiscI;
%....................................................................
% $$$ aveDiversityContI = aveEindexChapContI;
% $$$ aveDiversityDiscI = aveEindexChapDiscI;
%....................................................................
%====================================================================
% $$$ %....................................................................
% $$$ fignum = 1500;
% $$$ %....................................................................
% $$$ jamstecrest_loopruns_imagescsubplots_4x4(aveMUPcont,aveMUZcont,aveFPHYTcont,aveFZOOcont,avePHYTcont,aveZOOcont,aveDINcont,avePONcont,aveESDphyAveCont,aveDiversityCont,fignum,mypackages)
% $$$ %....................................................................
% $$$ fignum = 2500;
% $$$ %....................................................................
% $$$ jamstecrest_loopruns_imagescsubplots_4x4(aveMUPdisc,aveMUZdisc,aveFPHYTdisc,aveFZOOdisc,avePHYTdisc,aveZOOdisc,aveDINdisc,avePONdisc,aveESDphyAveDisc,aveDiversityDisc,fignum,mypackages)
% $$$ %....................................................................
%====================================================================
%....................................................................
fignum = 1510;
%....................................................................
[hfig] = jamstecrest_loopruns_imagescsubplots_4x4(aveMUPcontI,aveMUZcontI,aveFPHYTcontI,aveFZOOcontI,avePHYTcontI,aveZOOcontI,aveDINcontI,avePONcontI,aveESDphyAveContI,aveDiversityContI,fignum,mypackages)
%....................................................................
% $$$ print(hfig,'-dpng','-r200','jamstecrest_loopruns_imagescsubplots_4x4_KTW_press_cont.png')
% $$$ print(hfig,'-dpdf','-r200','jamstecrest_loopruns_imagescsubplots_4x4_KTW_press_cont.pdf')
% $$$ print(hfig,'-depsc','-r300','jamstecrest_loopruns_imagescsubplots_4x4_KTW_press_cont.eps')
%....................................................................
print(hfig,'-dpng','-r200','jamstecrest_loopruns_imagescsubplots_4x4_KTW_pulse_cont.png')
print(hfig,'-dpdf','-r200','jamstecrest_loopruns_imagescsubplots_4x4_KTW_pulse_cont.pdf')
print(hfig,'-depsc','-r300','jamstecrest_loopruns_imagescsubplots_4x4_KTW_pulse_cont.eps')
%....................................................................
fignum = 2510;
%....................................................................
jamstecrest_loopruns_imagescsubplots_4x4(aveMUPdiscI,aveMUZdiscI,aveFPHYTdiscI,aveFZOOdiscI,avePHYTdiscI,aveZOOdiscI,aveDINdiscI,avePONdiscI,aveESDphyAveDiscI,aveDiversityDiscI,fignum,mypackages)
%....................................................................
print('-dpng ','-r200','jamstecrest_loopruns_imagescsubplots_4x4_KTW_disc.png')
print('-depsc','-r300','jamstecrest_loopruns_imagescsubplots_4x4_KTW_disc.eps')
%....................................................................
%====================================================================
%....................................................................
fignum = 1000;
%....................................................................
myTitle001 = 'MUP';
myTitle002 = 'MUZ';
myTitle003 = 'SDIN';
myTitle004 = 'NTOT';
%....................................................................
jamstecrest_loopruns_imagescsubplots_2x2(aveMUPcontI,aveMUZcontI,aveSDINI,aveNTOTcontI,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages);
%....................................................................
fignum = 2000;
%....................................................................
myTitle001 = 'MUP';
myTitle002 = 'MUZ';
myTitle003 = 'SDIN';
myTitle004 = 'NTOT';
%....................................................................
jamstecrest_loopruns_imagescsubplots_2x2(aveMUPdiscI,aveMUZdiscI,aveSDINI,aveNTOTdiscI,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages);
%....................................................................
%====================================================================
%....................................................................
fignum = 1005;
%....................................................................
myTitle001 = 'MUP';
myTitle002 = 'MUZ';
myTitle003 = 'FPHYT';
myTitle004 = 'FZOO';
%....................................................................
jamstecrest_loopruns_imagescsubplots_2x2(aveMUPcontI,aveMUZcontI,aveFPHYTcontI,aveFZOOcontI,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages);
%....................................................................
fignum = 2005;
%....................................................................
myTitle001 = 'MUP';
myTitle002 = 'MUZ';
myTitle003 = 'FPHYT';
myTitle004 = 'FZOO';
%....................................................................
jamstecrest_loopruns_imagescsubplots_2x2(aveMUPdiscI,aveMUZdiscI,aveFPHYTdiscI,aveFZOOdiscI,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages);
%....................................................................
%====================================================================
%....................................................................
fignum = 1010;
%....................................................................
myTitle001 = 'PHY';
myTitle002 = 'ZOO';
myTitle003 = 'DIN';
myTitle004 = 'PON';
%....................................................................
jamstecrest_loopruns_imagescsubplots_2x2(avePHYTcontI,aveZOOcontI,aveDINcontI,avePONcontI,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages);
%....................................................................
fignum = 2010;
%....................................................................
myTitle001 = 'PHY';
myTitle002 = 'ZOO';
myTitle003 = 'DIN';
myTitle004 = 'PON';
%....................................................................
jamstecrest_loopruns_imagescsubplots_2x2(avePHYTdiscI,aveZOOdiscI,aveDINdiscI,avePONdiscI,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages);
%....................................................................
%====================================================================
%....................................................................
fignum = 1020;
%....................................................................
myTitle001 = 'logESD ave';
myTitle002 = 'logESD std';
myTitle003 = 'ESD ave';
myTitle004 = 'ESD cv';
%....................................................................
jamstecrest_loopruns_imagescsubplots_2x2(avelogESDphyAveContI,avelogESDphyStdContI,aveESDphyAveContI,aveESDphyCVcontI,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages);
%....................................................................
fignum = 2020;
%....................................................................
myTitle001 = 'logESD ave';
myTitle002 = 'logESD std';
myTitle003 = 'ESD ave';
myTitle004 = 'ESD cv';
%....................................................................
jamstecrest_loopruns_imagescsubplots_2x2(avelogESDphyAveDiscI,avelogESDphyStdDiscI,aveESDphyAveDiscI,aveESDphyCVdiscI,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages);
%....................................................................
%====================================================================
%....................................................................
fignum = 1030;
%....................................................................
myTitle001 = 'MUP';
myTitle002 = 'PHY';
myTitle003 = 'ESD ave';
myTitle004 = 'exp(H)';
%....................................................................
jamstecrest_loopruns_imagescsubplots_2x2(aveMUPcontI,avePHYTcontI,aveESDphyAveContI,aveEindexContI,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages);
%....................................................................
fignum = 2030;
%....................................................................
myTitle001 = 'MUP';
myTitle002 = 'PHY';
myTitle003 = 'ESD ave';
myTitle004 = 'exp(H)';
%....................................................................
jamstecrest_loopruns_imagescsubplots_2x2(aveMUPdiscI,avePHYTdiscI,aveESDphyAveDiscI,aveEindexDiscI,myTitle001,myTitle002,myTitle003,myTitle004,fignum,mypackages);
%....................................................................
%====================================================================
%....................................................................
%....................................................................
%....................................................................
%====================================================================

%********************************************************************
return


