function [Xrandn]=myRandn(nrows,ncols,xmean,sigma,seed1,seed2)

%********************************************************************
%Use: [x]=myRandn(nrows,ncols,mean,sigma,seed1,seed2)
%
%This function generates a Normal (ie. Gaussian) distribution 
%with values between 0.5 and 1.5 centered at 1.0 mean. 
%It is based on Box-Muller transform. 
%********************************************************************
%--------------------------------------------------------------------
%<http://answers.yahoo.com/question/index?qid=20080303124312AAMaYeP>
%Use the Box-Muller transform. 
%--------------------------------------------------------------------
%...........................................    
% $$$ xmean = 1; %enter the mean you want or need.
% $$$ sigma = 1/4; %enter the standard deviation you want or need.
%...........................................    
% $$$ xmean = 1; %enter the mean you want or need.
% $$$ sigma = 1/8; %enter the standard deviation you want or need.
%...........................................    
% $$$ xmean = 1; %enter the mean you want or need.
% $$$ sigma = 1/16; %enter the standard deviation you want or need.
%...........................................    
rand('seed',seed1);
u1 = rand(nrows,ncols);
%...........................................    
rand('seed',seed2);
u2 = rand(nrows,ncols);
%...........................................    
z1 = sqrt(-2.*log(u1)).*(sin(2*pi*u2));
z2 = sqrt(-2.*log(u1)).*(cos(2*pi*u2));
%...........................................    
x1 = xmean + z1*sigma;
x2 = xmean + z2*sigma;
%...........................................    
% $$$ x = x1; %original.
%...........................................    
x = (x1+x2)/2; %mio (no se si esta del todo bien).
%...........................................
Hinf = find(x<=0);
Hsup = find(x>=2);
%...........................................
x(Hinf) = xmean;
x(Hsup) = xmean;
%...........................................
%%%%%%%%
%OUTPUT:
%%%%%%%%
Xrandn = x;

