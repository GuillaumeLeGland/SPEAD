function [ShannonEntropyTheoreticalAbsnormal1D,ShannonEntropyTheoreticalLognormal1D] = jamstecrest_ShannonContinuousEntropy1D(Xave,Xstd,fignum,mypackages)

%========================================================================
%MY PACKAGES FOR PLOTING:
%........................................................................
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal; 
%........................................................................
%========================================================================
%........................................................................
ndim = 1;
%........................................................................
xaverg = Xave;
xsigma = Xstd;
%........................................................................
ShannonEntropyTheoreticalAbsnormal1D = (ndim/2) * (1 + log(2*pi)) + (1/2)*log(xsigma.^2);
%........................................................................
% <https://en.wikipedia.org/wiki/Log-normal_distribution> 
ShannonEntropyTheoreticalLognormal1D = log(xsigma .* exp(xaverg + 1/2) * sqrt(2*pi));
%........................................................................
% Quintana etal 2008:
ShannonEntropyTheoreticalAbsnormal1Dbis = (1/2) + log(sqrt(2*pi)*xsigma);
ShannonEntropyTheoreticalLognormal1Dbis = (1/2) + log(sqrt(2*pi)*xsigma) + xaverg;
%........................................................................
figure(fignum)
subplot(2,2,1)
disp(size(ShannonEntropyTheoreticalAbsnormal1D))
imagesc(ShannonEntropyTheoreticalAbsnormal1D)
title('Entropy of Abs Normal Distribution')
colorbar_funhan(verticales)
subplot(2,2,2)
imagesc(ShannonEntropyTheoreticalLognormal1D)
title('Entropy of Log Normal Distribution')
colorbar_funhan(verticales)
grid on
subplot(2,2,3)
imagesc(ShannonEntropyTheoreticalAbsnormal1Dbis)
title('Entropy of Abs Normal Distribution')
colorbar_funhan(verticales)
grid on
subplot(2,2,4)
imagesc(ShannonEntropyTheoreticalLognormal1Dbis)
title('Entropy of Log Normal Distribution')
colorbar_funhan(verticales)
grid on
%........................................................................
%========================================================================
%************************************************************************
return
