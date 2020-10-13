function [XweightAve,XweightStd,YweightAve,YweightStd,XYweightCor] = SPEAD_1D_covariance(PHY,XTRAIT,YTRAIT,nxphy,nyphy)
% Function to plot weighted means, standard deviations and correlations of the traits

%===================================================================================
% <http://en.wikipedia.org/wiki/Weighted_arithmetic_mean> 
% <http://en.wikipedia.org/wiki/Mean_square_weighted_deviation> 
%...................................................................................
sumPHY = sum(PHY(:,:,:),3); 
sumPHYrepmat = repmat(sumPHY,[1,1,nxphy,nyphy]);
Weights = PHY ./ sumPHYrepmat;
sumWeights = sum(Weights(:,:,:),3);
%...................................................................................
XTR = repmat(XTRAIT,[1,1,1,nyphy]);
YTR = permute(repmat(YTRAIT,[1,1,1,nxphy]),[1,2,4,3]);
%...................................................................................
Xweight = Weights .* XTR;
Yweight = Weights .* YTR;
%...................................................................................
XweightSum = sum(Xweight(:,:,:),3);
YweightSum = sum(Yweight(:,:,:),3);
XweightAve = XweightSum ./ sumWeights;
YweightAve = YweightSum ./ sumWeights;
%...................................................................................
XweightAveRepmat = repmat(XweightAve,[1,1,nxphy,nyphy]);
Xdispersion = Weights .* (XTR - XweightAveRepmat).^2;
sumXdispersion = sum(Xdispersion(:,:,:),3);
XweightVar = sumXdispersion ./ sumWeights;
XweightStd = sqrt(XweightVar);
%...................................................................................
YweightAveRepmat = repmat(YweightAve,[1,1,nxphy,nyphy]);
Ydispersion = Weights .* (YTR - YweightAveRepmat).^2;
sumYdispersion = sum(Ydispersion(:,:,:),3);
YweightVar = sumYdispersion ./ sumWeights;
YweightStd = sqrt(YweightVar);
%...................................................................................
XYdispersion = Weights .* (XTR - XweightAveRepmat) .* (YTR - YweightAveRepmat);
sumXYdispersion = sum(XYdispersion(:,:,:),3);
XYweightCov = sumXYdispersion ./ sumWeights;
XYweightCor = XYweightCov ./ (XweightStd .* YweightStd);
%...................................................................................
%***********************************************************************************
return




