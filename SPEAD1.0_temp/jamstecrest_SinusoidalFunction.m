function [fsin] = jamstecrest_SinusoidalFunction(ywin,ysum,ndays,tmax,shift,keySinusoidalShape)
%function [fsin] = jamstecrest_SinusoidalFunction(ywin,ysum,ndays,tmax,keySinusoidalShape)
%function [fsin] = jamstecrest_SinusoidalFunction(ymin,ymax,ndays,tmax,keySinusoidalShape)

%........................................................................
T = ndays;
tau = ndays/4;
time = [1:1:tmax];
%ttime = time - tau;
ttime = time - tau - shift; % Shift compared with PAR (Le Gland, 27/09/2019)
wrad = (2*pi) / T;
%........................................................................
%A0 = ymin;
%A = (ymax-ymin)/2; %amplitud.
A0 = ywin; % Extrema are specified in term of "winter" and "summer" (Le Gland, 26/07/2019)
A = (ysum-ywin)/2;
ysin = A0 + A*(1.0 + sin(wrad*ttime)); %Sinusoidal function.
%........................................................................
if strcmp(keySinusoidalShape,'Linear')
    fsin = ysin; 
elseif strcmp(keySinusoidalShape,'Quadratic')
    npower = 1.5;
    % xsin = (ysin/ymax).^2 / (ymin./ymax).^2; % ymax is the min and ymin is the max ! (Le Gland, 26/07/2019)
    xsin = ((ysin/ysum).^2 - 1) ./ ((ywin./ysum).^2 - 1); % In order to stay between 0 and 1 (otherwise the minimum is never exactly reached)
    fsin = ysum + 2*(-A)*xsin.^npower;  
elseif strcmp(keySinusoidalShape,'Step')
    % Step function, representing violent perturbation (Le Gland, 29/04/2019)
    fsin = ymin + A*(1.0 + sign(sin(26*wrad*ttime)));
end
%........................................................................
return
