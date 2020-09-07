function [fsin] = mySinusoidalFunction(ymin,ymax,ndays,tmax,tau,keySinusoidalShape)

%........................................................................
T = ndays;
%%tau = ndays/4;
time = [1:1:tmax];
ttime = time - tau;
wrad = (2*pi) / T;
%........................................................................
A0 = ymin;
A = (ymax-ymin)/2; %amplitud.
ysin = A0 + A*(1.0 + sin(wrad*ttime)); %Sinusoidal function.
%........................................................................
if strcmp(keySinusoidalShape,'Linear')
    fsin = ysin; 
elseif strcmp(keySinusoidalShape,'Quadratic')
    npower = 1.5;
    xsin = (ysin/ymax).^2 / (ymin./ymax).^2;
    fsin = ymax + 2*(-A)*xsin.^npower;
end
%........................................................................

