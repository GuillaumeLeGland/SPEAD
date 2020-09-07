function [fx,xaxis,xmean,xsigma] = myGaussianDistribution1D(xmin,xmax,sigmaPcnt,nbins)

%===================================================================================
%...................................................................................
xdel = ((xmax-xmin)/(nbins-1)); 
xaxis = [xmin:xdel:xmax]; %[log(um)] 
xmean = mean([xmin:xmax]);
xsigma = sigmaPcnt*(xmax-xmin); %Okay. 
fx = (1.0 / (xsigma * sqrt(2*pi))) * exp( -(xaxis - xmean).^2 / (2*xsigma^2) ); %Okay.
%...................................................................................
%===================================================================================

return
