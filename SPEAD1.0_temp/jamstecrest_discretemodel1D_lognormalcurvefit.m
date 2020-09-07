function [logESDave,logESDstd,B1coeff,B2coeff,R2coeff,ShannonEntropy,ShannonEntropyTheoretical] = jamstecrest_discretemodel1D_lognormalcurvefit(PHYsspdisc3D,ESDphy,logESDphy,xsizedel)

%========================================================================
disp(' jamstecrest_discretemodel1D_lognormalcurvefit.m -- START **')
%........................................................................
maxPHY = max(PHYsspdisc3D(:)); %[depth,time,species]
%........................................................................
% [mdepths,mdays,mspecies] = size(PHYsspdisc3D); 
[mdepths,mdays,mspecies] = size(PHYsspdisc3D(:,:,:)); % 2 traits (Le Gland, 18/07/2019)
%........................................................................
%========================================================================
%........................................................................
for zdepthj = 1:mdepths 
    zdepthj 
    for dayj = 1:mdays
	%................................................................
	%phy = squeeze(PHYsspdisc3D(zdepthj,dayj,:));
  % Case with 2 traits (Le Gland, 18/07/2019)
  phy2D = squeeze(PHYsspdisc3D(zdepthj,dayj,:,:));
  phy = sum(phy2D,2); % sum over all values of optimal temperature
	phytot = sum(phy);
	phypcnt = phy./phytot;
	%................................................................
	%================================================================
	%................................................................
% $$$ 	ESDweighted = phypcnt .* ESDphy(:); 
% $$$ 	ESDgeomean = sum(ESDweighted); %Geometric mean. 
% $$$ 	ESDdispersion = phypcnt .* (ESDphy(:) - ESDgeomean).^2; 
% $$$ 	%................................................................
% $$$ 	ESDvariance = sum(ESDdispersion); 
% $$$ 	ESDsigma = sqrt(ESDvariance); 
	%................................................................
	%================================================================
	%NOTE: I THINK I SHOULD BE CALLING THIS "WEIGTHED ARITMETIC MEAN" 
	%INSTEAD OF "GEOMETRIC MEAN" WHICH IS ANOTHER DIFFERENT THING!!!
	%----------------------------------------------------------------
	% <https://en.wikipedia.org/wiki/Weighted_arithmetic_mean> 
	%----------------------------------------------------------------
	%................................................................
	logESDweighted = phypcnt .* logESDphy(:); 
	logESDgeomean = sum(logESDweighted); %Geometric mean. 
	logESDdispersion = phypcnt .* (logESDphy(:) - logESDgeomean).^2; 
	%................................................................
	logESDvariance = sum(logESDdispersion); 
	logESDsigma = sqrt(logESDvariance); 
	%................................................................
	%================================================================
	%MATLAB FUNCTIONS:
	%................................................................
% $$$ 	logESDgeomeanBis = wmean(logESDphy(:),phypcnt(:));
% $$$ 	%................................................................
% $$$ 	logESDsigmaBis = std(logESDphy(:),phypcnt(:));
	%................................................................
	%================================================================
	%................................................................
	xdat = logESDphy; 
	xave = logESDgeomean; 
	xsig = logESDsigma; 
	%................................................................
	PDF = (1.0 / (xsig * sqrt(2*pi))) * exp( -(xdat - xave).^2 / (2*xsig^2) ); %Okay.
	%................................................................
	sumPDF = sum(PDF .* xsizedel); %must be equal to one!!! 
	%................................................................
	phychap = phytot .* (PDF .* xsizedel); 
	%................................................................
	%================================================================
	%LINEAR REGRESSION: 
	%................................................................
	x = phy(:); 
	y = phychap(:); 
	%................................................................
	xr = x; 
	%................................................................
	X = [ones(length(x),1),x];
	Xr =[ones(length(xr),1),xr];
	%................................................................
	B = (X'*X)\(X'*y); %coefficientes de regression.
	%................................................................
	ychap = Xr*B; %recta linear de regression.
	%................................................................
	ym = sum(y)/length(y);
	%................................................................
	vt  = sum((y-ym).^2); 
	ve  = sum((ychap-ym).^2); 
	vne = sum((y-ychap).^2); 
	%................................................................
	r2 = ve/vt; %Coefficient of determination. 
	%................................................................
	%================================================================
	%................................................................
% $$$ 	figure(10)
% $$$ 	dayj
% $$$ 	subplot(2,2,1) 
% $$$ 	plot(logESDphy,phy,'r*')
% $$$ 	hold on
% $$$ 	plot(logESDphy,phychap,'-b.')
% $$$ 	hold off
% $$$ 	grid on
% $$$ 	subplot(2,2,2)
% $$$ 	plot(phy,phychap,'*')
% $$$ 	hold on
% $$$ 	plot(x,xr)
% $$$ 	hold off 
% $$$ 	axis([0 maxPHY, 0 maxPHY]) 
% $$$ 	grid on
% $$$ 	pause(0.1)
	%................................................................
	%================================================================
	%................................................................
	logESDave(zdepthj,dayj) = logESDgeomean; 
	logESDstd(zdepthj,dayj) = logESDsigma; 
	%................................................................
	B1coeff(zdepthj,dayj) = B(1); 
	B2coeff(zdepthj,dayj) = B(2); 
	R2coeff(zdepthj,dayj) = r2; 
	%................................................................
	%================================================================
	%................................................................
	%%[ijShannonEntropy,ijShannonEntropyTheoretical] = myShannonContinousEntropyFromLanSmith(phy,logESDphy,logESDgeomean,logESDsigma);
	%................................................................
	[ijShannonEntropy001,ijShannonEntropyTheoretical001] = myShannonContinousEntropyFromLanSmith(phychap,logESDphy,logESDgeomean,logESDsigma);
	%................................................................
	[ijShannonEntropy002,ijShannonEntropyTheoretical002] = myShannonContinousEntropy1D(phychap,logESDphy,logESDgeomean,logESDsigma);
	%................................................................
	%================================================================
	%................................................................
% $$$ 	ShannonEntropy(zdepthj,dayj) = ijShannonEntropy001;
% $$$ 	ShannonEntropyTheoretical(zdepthj,dayj) = ijShannonEntropyTheoretical001;
	%................................................................
	ShannonEntropy(zdepthj,dayj) = ijShannonEntropy002;
	ShannonEntropyTheoretical(zdepthj,dayj) = ijShannonEntropyTheoretical002;
	%................................................................
	%================================================================
    end
end
%........................................................................
%========================================================================
disp(' jamstecrest_discretemodel1D_lognormalcurvefit.m -- END **')
%***********************************************************
return 

%========================================================================
%........................................................................
xmean = 0;
xsigma = 1;
%........................................................................
ypdf = lognpdf(ESDphy,xmean,xsigma); 
%........................................................................
figure(100)
subplot(2,2,1) 
plot(logESDphy,ypdf,'-')
hold on
plot(logESDphy,ypdf,'r*')
hold off
grid on
%........................................................................
%========================================================================
