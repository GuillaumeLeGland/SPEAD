function [ShannonEntropyNumerically,ShannonEntropyTheoretical] = myShannonContinousEntrophy(phy,logESD,logESDave,logESDstd);
% <https://www.hiit.fi/u/ahonkela/dippa/node94.html>
% <https://en.wikipedia.org/wiki/Differential_entropy> 
% <http://stackoverflow.com/questions/16527672/entropy-estimation-using-histogram-of-normal-data-vs-direct-formula-matlab> 

%===================================================================================
%................................................................................... 
p = phy./nansum(phy); 
dx = diff(logESD); 
dx = [dx,dx(end)]; %add extra point to match size.
%................................................................................... 
ShannonEntropyClassical = -nansum( p.*log(p)); 
%................................................................................... 
ShannonEntropyHistogram = -nansum((p.*log(p)).*dx); 
%................................................................................... 
ShannonEntropyTheoreticalAbsnormal002 = log(sqrt(2*pi*exp(1))*logESDstd); %Should be the same.
ShannonEntropyTheoreticalAbsnormal001 = (1/2)*log(2*pi*exp(1)*logESDstd.^2);
ShannonEntropyTheoreticalAbsnormal003 = (1/2) * log((2*pi)*(logESDstd.^2)) + (1/2);
%................................................................................... 
ShannonEntropyTheoreticalLognormal001 = (1/2) + log(sqrt(2*pi) * logESDstd) + logESDave;  %Quintana.
ShannonEntropyTheoreticalLognormal002 = (1/2) * log(2*pi*exp(1)* logESDstd^2) + logESDave; %Wikipedia.
%................................................................................... 
%===================================================================================
%................................................................................... 
nbins = length(phy); 
%................................................................................... 
peven = (1/nbins) * ones(1,nbins); 
%................................................................................... 
ShannonEntropyClassicalMax = -nansum( (peven) .* log(peven) );
ShannonEntropyHistogramMax = -nansum( (peven) .* log(peven) .*dx );
%................................................................................... 
ShannonEntropyClassicalStar = ShannonEntropyClassical/ShannonEntropyClassicalMax;
ShannonEntropyHistogramStar = ShannonEntropyHistogram/ShannonEntropyHistogramMax;
%................................................................................... 
%===================================================================================
%................................................................................... 
delta = (max(logESD)-min(logESD))/(length(logESD)-1);
lower = min(logESD)-delta/2;
upper = max(logESD)+delta/2;
ncell = length(logESD);
descriptor = [lower,upper,ncell];
ShannonEntropyNumericallyBis = myentropy(phy,descriptor);
%................................................................................... 
%===================================================================================
%................................................................................... 
%%ShannonEntropyNumerically = ShannonEntropyClassical;
ShannonEntropyNumerically = ShannonEntropyHistogram;
%................................................................................... 
%%ShannonEntropyNumerically = ShannonEntropyClassicalStar;
%%ShannonEntropyNumerically = ShannonEntropyHistogramStar;
%................................................................................... 
ShannonEntropyTheoretical = ShannonEntropyTheoreticalAbsnormal003;
%%ShannonEntropyTheoretical = ShannonEntropyTheoreticalLognormal002;
%................................................................................... 
%===================================================================================

% $$$ ShannonEntropy
% $$$ ShannonEntropyBis 
% $$$ ShannonEntropyTheoretical

%**********************************************************************
return

%===================================================================================
%................................................................................... 
xnbins = 64;
ynbins = 64;
%................................................................................... 
xmin = -1; %For Logn (ie. ESDmin = 0.36788) 
ymin = -1; %For Logn (ie. ESDmin = 0.36788) 
%................................................................................... 
xmax = +4; %For Logn (ie. ESDmax = 54.598) 
ymax = +4; %For Logn (ie. ESDmax = 54.598) 
%................................................................................... 
xsigmaPcnt = 0.2; 
ysigmaPcnt = 0.1;
%................................................................................... 
[fx,xaxis,xmean,xsigma] = myGaussianDistribution1D(xmin,xmax,xsigmaPcnt,xnbins);
[fy,yaxis,ymean,ysigma] = myGaussianDistribution1D(ymin,ymax,ysigmaPcnt,ynbins);
%................................................................................... 
xdel = diff(xaxis);
ydel = diff(yaxis);
%................................................................................... 
xdel = [xdel,xdel(1)];
ydel = [ydel,ydel(1)];
%................................................................................... 
xlag = 0;
ylag = 0;
%................................................................................... 
%===================================================================================
%................................................................................... 
[fx64] = fx;
[fx32] = mychgresolvector(fx64,2);
[fx16] = mychgresolvector(fx32,2);
[fx08] = mychgresolvector(fx16,2);
%................................................................................... 
[xaxis64] = xaxis;
[xaxis32] = mychgresolvector(xaxis64,2);
[xaxis16] = mychgresolvector(xaxis32,2);
[xaxis08] = mychgresolvector(xaxis16,2);
%................................................................................... 
xdel64 = mean(diff(xaxis64));
xdel32 = mean(diff(xaxis32));
xdel16 = mean(diff(xaxis16));
xdel08 = mean(diff(xaxis08));
%................................................................................... 
sumfx64 = sum(fx64 .* xdel64); %must be equal to one!!! 
sumfx32 = sum(fx32 .* xdel32); %must be equal to one!!! 
sumfx16 = sum(fx16 .* xdel16); %must be equal to one!!! 
sumfx08 = sum(fx08 .* xdel08); %must be equal to one!!! 
%................................................................................... 
%===================================================================================
%................................................................................... 
[fy64] = fy;
[fy32] = mychgresolvector(fy64,2);
[fy16] = mychgresolvector(fy32,2);
[fy08] = mychgresolvector(fy16,2);
%................................................................................... 
[yaxis64] = yaxis;
[yaxis32] = mychgresolvector(yaxis64,2);
[yaxis16] = mychgresolvector(yaxis32,2);
[yaxis08] = mychgresolvector(yaxis16,2);
%................................................................................... 
ydel64 = mean(diff(yaxis64));
ydel32 = mean(diff(yaxis32));
ydel16 = mean(diff(yaxis16));
ydel08 = mean(diff(yaxis08));
%................................................................................... 
sumfy64 = sum(fy64 .* ydel64); %must be equal to one!!! 
sumfy32 = sum(fy32 .* ydel32); %must be equal to one!!! 
sumfy16 = sum(fy16 .* ydel16); %must be equal to one!!! 
sumfy08 = sum(fy08 .* ydel08); %must be equal to one!!! 
%................................................................................... 
%===================================================================================
%................................................................................... 
figure(10)
subplot(2,2,1)
plot(xaxis64,fx64,'b-',xaxis64,fx64,'r.')
subplot(2,2,2)
plot(xaxis32,fx32,'b-',xaxis32,fx32,'r.')
subplot(2,2,3)
plot(xaxis16,fx16,'b-',xaxis16,fx16,'r.')
subplot(2,2,4)
plot(xaxis08,fx08,'b-',xaxis08,fx08,'r.')
%................................................................................... 
figure(20)
subplot(2,2,1)
plot(yaxis64,fy64,'b-',yaxis64,fy64,'r.')
subplot(2,2,2)
plot(yaxis32,fy32,'b-',yaxis32,fy32,'r.')
subplot(2,2,3)
plot(yaxis16,fy16,'b-',yaxis16,fy16,'r.')
subplot(2,2,4)
plot(yaxis08,fy08,'b-',yaxis08,fy08,'r.')
%................................................................................... 
%===================================================================================
%................................................................................... 
[ShannonEntropy64,ShannonEntropyTheoretical64] = myShannonContinousEntrophy(fx64,xaxis64,xmean,xsigma);
[ShannonEntropy32,ShannonEntropyTheoretical32] = myShannonContinousEntrophy(fx32,xaxis32,xmean,xsigma);
[ShannonEntropy16,ShannonEntropyTheoretical16] = myShannonContinousEntrophy(fx16,xaxis16,xmean,xsigma);
[ShannonEntropy08,ShannonEntropyTheoretical08] = myShannonContinousEntrophy(fx08,xaxis08,xmean,xsigma);
%................................................................................... 
[ShannonEntropy64,ShannonEntropyTheoretical64]
[ShannonEntropy32,ShannonEntropyTheoretical32]
[ShannonEntropy16,ShannonEntropyTheoretical16]
[ShannonEntropy08,ShannonEntropyTheoretical08]
%................................................................................... 
%===================================================================================
%................................................................................... 
[ShannonEntropyND64fx,ShannonEntropyTheoreticalND64fx] = myShannonContinousEntrophyND(xaxis64(:),fx64(:));
[ShannonEntropyND32fx,ShannonEntropyTheoreticalND32fx] = myShannonContinousEntrophyND(xaxis32(:),fx32(:));
[ShannonEntropyND16fx,ShannonEntropyTheoreticalND16fx] = myShannonContinousEntrophyND(xaxis16(:),fx16(:));
[ShannonEntropyND08fx,ShannonEntropyTheoreticalND08fx] = myShannonContinousEntrophyND(xaxis08(:),fx08(:));
%................................................................................... 
[ShannonEntropyND64fx,ShannonEntropyTheoreticalND64fx]
[ShannonEntropyND32fx,ShannonEntropyTheoreticalND32fx]
[ShannonEntropyND16fx,ShannonEntropyTheoreticalND16fx]
[ShannonEntropyND08fx,ShannonEntropyTheoreticalND08fx]
%................................................................................... 
%===================================================================================
%................................................................................... 
[ShannonEntropyND64fxx,ShannonEntropyTheoreticalND64fxx] = myShannonContinousEntrophyND(xaxis64(:),fx64(:),xaxis64(:),fx64(:));
[ShannonEntropyND32fxx,ShannonEntropyTheoreticalND32fxx] = myShannonContinousEntrophyND(xaxis32(:),fx32(:),xaxis32(:),fx32(:));
[ShannonEntropyND16fxx,ShannonEntropyTheoreticalND16fxx] = myShannonContinousEntrophyND(xaxis16(:),fx16(:),xaxis16(:),fx16(:));
[ShannonEntropyND08fxx,ShannonEntropyTheoreticalND08fxx] = myShannonContinousEntrophyND(xaxis08(:),fx08(:),xaxis08(:),fx08(:));
%................................................................................... 
[ShannonEntropyND64fxx,ShannonEntropyTheoreticalND64fxx]
[ShannonEntropyND32fxx,ShannonEntropyTheoreticalND32fxx]
[ShannonEntropyND16fxx,ShannonEntropyTheoreticalND16fxx]
[ShannonEntropyND08fxx,ShannonEntropyTheoreticalND08fxx]
%................................................................................... 
%===================================================================================
%................................................................................... 
[ShannonEntropyND64fxy,ShannonEntropyTheoreticalND64fxy] = myShannonContinousEntrophyND(xaxis64(:),fx64(:),yaxis64(:),fy64(:));
[ShannonEntropyND32fxy,ShannonEntropyTheoreticalND32fxy] = myShannonContinousEntrophyND(xaxis32(:),fx32(:),yaxis32(:),fy32(:));
[ShannonEntropyND16fxy,ShannonEntropyTheoreticalND16fxy] = myShannonContinousEntrophyND(xaxis16(:),fx16(:),yaxis16(:),fy16(:));
[ShannonEntropyND08fxy,ShannonEntropyTheoreticalND08fxy] = myShannonContinousEntrophyND(xaxis08(:),fx08(:),yaxis08(:),fy08(:));
%................................................................................... 
[ShannonEntropyND64fxy,ShannonEntropyTheoreticalND64fxy]
[ShannonEntropyND32fxy,ShannonEntropyTheoreticalND32fxy]
[ShannonEntropyND16fxy,ShannonEntropyTheoreticalND16fxy]
[ShannonEntropyND08fxy,ShannonEntropyTheoreticalND08fxy]
%................................................................................... 
%===================================================================================
%................................................................................... 
Ratiofxx_fx  = ShannonEntropyTheoreticalND64fxy ./ ShannonEntropyTheoreticalND64fx
Ratiofxy_fx  = ShannonEntropyTheoreticalND64fxy ./ ShannonEntropyTheoreticalND64fx
Ratiofxy_fxx = ShannonEntropyTheoreticalND64fxy ./ ShannonEntropyTheoreticalND64fxx
%................................................................................... 
%===================================================================================
%................................................................................... 
xsigma = xsigmaPcnt*(xmax-xmin); %Okay. 
ysigma = ysigmaPcnt*(ymax-ymin); %Okay. 
%................................................................................... 
[FXX,XAXIS,YAXIS] = myGaussianDistribution2D(xmin,xmin,xmax,xmax,xlag,xlag,xsigma,xsigma,xnbins,xnbins);
%................................................................................... 
[FXY,XAXIS,YAXIS] = myGaussianDistribution2D(xmin,ymin,xmax,ymax,xlag,ylag,xsigma,ysigma,xnbins,ynbins);
%................................................................................... 
FXY001 = 1/(2*pi*sigmax*sigmay) * exp(-( (XAXIS - xm).^2 / (2*sigmax^2) + (YAXIS - ym).^2 / (2*sigmay^2) ));
%................................................................................... 
%===================================================================================


nbins = 64; 
samples = 1000; 
[p,x] = hist(samples,nbins);
area = (x(2)-x(1))*sum(p);
p = p/area;

dp = (x(2)-x(1));
area = sum(p)*dp;
H_histogram = (1/2)*log2(2*pi*exp(1))
H_theoretical = -nansum((p*dp).*log2(p))

