function [] = myGaussianAsymmetric1D()

xaxis = [0:0.1:36];
xmean = 16;
xsigma = 4;
factor = 1.1;
%%factor = 1.5;
lambda = factor*xmean; 
sigmafilter = 1.0;
sigmafilter = 1.0 - xaxis / lambda; %Asymmetric Gaussian. %OKAY 

fx_sym001  = exp(-(xaxis - xmean).^2 ./  (2*xsigma)^2);
fx_asym001 = exp(-(xaxis - xmean).^2 ./ ((2*xsigma*sigmafilter).^2));

fx_sym002  = exp(-(xaxis - xmean).^2 ./ (2*xsigma^2));
fx_asym002 = exp(-(xaxis - xmean).^2 ./ (2*xsigma^2 * sigmafilter.^2));

fx_sym003  = exp(-(xaxis - xmean).^2 ./ (2*xsigma^2));
fx_asym003 = exp(-(xaxis - xmean).^2 ./ (2*xsigma^2 * 2*sigmafilter.^2));
%%fx_asym003 = exp(-(xaxis - xmean).^2 ./ (2*xsigma^2 * (2*sigmafilter).^2));

figure(1)
subplot(2,2,1)
plot(xaxis,fx_sym001,'.b-')
hold on
plot(xaxis,fx_asym001,'.r-')
hold off
grid on
subplot(2,2,2)
plot(xaxis,fx_sym002,'.b-')
hold on
plot(xaxis,fx_asym002,'.r-')
hold off
grid on
subplot(2,2,3)
plot(xaxis,fx_sym003,'.b-')
hold on
plot(xaxis,fx_asym003,'.r-')
hold off
grid on

return

