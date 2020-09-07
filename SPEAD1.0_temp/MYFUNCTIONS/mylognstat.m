function [XmedianUpper,XmedianLower,Xmedian,Xstd,Xcvd,Xmean,Xvar,Xsig] = mylognstat(LOGXave,LOGXstd)

%===================================================================================
%-----------------------------------------------------------------------------------
%NOTE: 
%Doing ESDave = exp(logESDave) is okay. 
%Doing ESDstd = exp(logESDstd) is wrong!!!! (see how to do it right below) 
%-----------------------------------------------------------------------------------
%===================================================================================
%...................................................................................
%COMPUTING "MEDIAN" AND "RANGE" OF ABSOLUTE DISTRIBUTION FROM LOGNORMAL DISTRIBUTION "MEAN" AND "SIGMA":
Xmedian = exp(LOGXave); %Okay. 
XmedianUpper = Xmedian .* exp(LOGXstd);
XmedianLower = Xmedian ./ exp(LOGXstd);
%...................................................................................
%===================================================================================
%COMPUTING COEFFICIENT OF VARIATION OF CELL SIZE IN ABSOLUTE SCALE FROM LOG(CELL SIZE) SCALE:
%-----------------------------------------------------------------------------------
%NOTE: Xstd is not really the "standard deviation" but half the "delta"
%between the maximum and the minimum value of x in absolute distribution 
%of the same area covered by 2*logxsigma in the lognormal distribution. 
%This derivation comes from Lan Smith (JAMSTEC). See my blue libreta. 
%-----------------------------------------------------------------------------------
%...................................................................................
Xdelta = (XmedianUpper - XmedianLower); %Crude estimation but should be very accurate.
%...................................................................................
XdeltaChap001 = (exp(LOGXave) ./ exp(LOGXstd)) .* (exp(LOGXstd).^2 - 1d0); %Theoretical equation. 
XdeltaChap002 = (exp(LOGXave) ./ exp(LOGXstd)) .* (exp(2*LOGXstd)  - 1d0); %Theoretical equation. 
XdeltaChap = XdeltaChap001; 
%...................................................................................
%%Xstd = Xdelta; %STD (not sure if is better to use this or to use Delta / 2)
Xstd = Xdelta / 2; %(maybe STD is half the Delta ??)
%...................................................................................
Xcvd = Xstd./Xmedian; %Coefficient of variation (std/mean)
%...................................................................................
%===================================================================================
%COMPUTE "MEAN" AND "VARIANCE" OF ABSOLUTE DISTRIBUTION FROM LOGNORMAL DISTRIBUTION "MEAN" AND "SIGMA":
%...................................................................................
logx_ave = LOGXave;
%%logx_sig = LOGXstd; %Original. (okay when using above: Xstd = Xdelta / 2)
logx_sig = sqrt(LOGXstd); %Sale mas razonable pero no se por que? (tal vez porque estaba usando Xstd = Xdelta en lugar the Xstd = Xdelta/2 ???)
%...................................................................................
Xmean = exp(  logx_ave + (logx_sig.^2)/2);
Xvar =  exp(2*logx_ave + (logx_sig.^2)) .* (exp(logx_sig.^2) - 1d0);
Xsig = sqrt(Xvar); 
%...................................................................................
%===================================================================================
%...................................................................................
% $$$ Xmedian 
% $$$ Xcvd 
% $$$ Xstd 
% $$$ XstdChap 
%...................................................................................
% $$$ Xmean 
% $$$ Xvar 
% $$$ Xsig 
%...................................................................................
%===================================================================================
return
