function [] = jamstecrest_gaussecomodel1D_lognormal(ESD) %DOES NOT WORK!!!!!!
% <https://statisticaleconomics.files.wordpress.com/2014/01/various-properties-of-the-lognormal-distribution.pdf>

x = ESD;
logxave = 1.5;
logxsig = 0.5;

[xave,xsig] = lognstat(logxave,logxsig);

xvar = xsig.^2;

logxvar = log(1 + (xvar/xave^2));
logxsigBis = sqrt(logxvar)

xsigBis = sqrt(log(1 + xvar/xave^2));
xpsi = log(xave) - (1/2)*xsig^2;

fx = (1.0 ./ (x * sqrt(2*pi*xsig^2))) .* exp( (log(x) - xpsi).^2 / (2*xsig^2));

fx = (1.0 / (xsigma * sqrt(2*pi))) * exp( -(xtrait - xmean).^2 / (2*xsigma^2) );


H = (1/2) * log(2*pi*exp(1)*xsig^2) + xpsi
return
%***********************************************
