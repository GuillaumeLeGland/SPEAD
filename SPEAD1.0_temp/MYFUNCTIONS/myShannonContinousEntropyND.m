function [ShannonEntropyAbsnormalND,ShannonEntropyTheoreticalAbsnormalND] = myShannonContinousEntrophyND(varargin)

%***************************************************************************
%Use: function [] = myShannonContinousEntrophyND(xaxis,fx,yaxis,fy,zaxis,fz,...,naxis,fn)
%
%---------------------------------------------------------------------------
% <https://en.wikipedia.org/wiki/Multivariate_normal_distribution>
% 
% Entropy absnormal ND:
%---------------------------------------------------------------------------
% H = (ndim/2) * (1 + log(2*pi)) + (1/2)*log(det(xcov)) %OKAY!!!
%---------------------------------------------------------------------------
% H = (ndim/2) * log((2*pi)*exp(1)*(det(xcov))).^(1/ndim) %WRONG
%---------------------------------------------------------------------------
% H = (ndim/2) * ln((2*pi)*exp(1)*(sigmax^2 * sigmay^2 * ... * sigman^2))^(1/ndim)  %WRONG
%---------------------------------------------------------------------------
%***************************************************************************
%===========================================================================
%...........................................................................
groups = 2; %pairs of data in varargin [axis,fx]
nlen = length(varargin);
ndim = nlen/groups; 
xcov = zeros(ndim,ndim);
%...........................................................................
j = 0;
for index = [1:groups:nlen] 
    j = j + 1; 
    %.......................................................................
    xaxis = varargin{index}; 
    fx    = varargin{index+1}; 
    %.......................................................................
    px = fx / sum(fx);
    dx = diff(xaxis);
    dx = [dx;dx(end)]; %add extra point to match size.
    %.......................................................................
    xavej = sum(xaxis.*px); %mean of normal distribution-j
    xvarj = sum(px.*(xaxis-xavej).^2); %variance of normal distribution-j
    xstdj = sqrt(xvarj); %standard deviation of normal distribution-j
    %.......................................................................
    FX(:,j)   = fx(:);  %Array made of column vectors of all normal distributions-j
    xave(j)   = xavej;  %Vector of elements with the MEAN of all normal distributions-j
    xcov(j,j) = xvarj;  %Array of diagonal elements with the VARIANCE of all normal distributions-j 
    %.......................................................................
% $$$     jShannonEntropyTheoreticalAbsnormal1D_001 = (1/2) * log((2*pi)*(xstdj.^2)) + (1/2);
% $$$     jShannonEntropyTheoreticalAbsnormal1D_002 = (1/2) * log(2*pi*exp(1)*xstdj^2);
% $$$     jShannonEntropyHistogram = -nansum((px.*log(px)).*dx);
    %.......................................................................
end
%...........................................................................
ShannonEntropyAbsnormalND = nan; %DONT KNOW HOW TO COMPUTE THAT WITH CLASSICAL EQUATION BUT IN 2D.
%...........................................................................
ShannonEntropyTheoreticalAbsnormalND = (ndim/2) * (1 + log(2*pi)) + (1/2)*log(det(xcov))
%...........................................................................
%===========================================================================
%***************************************************************************
return
