function [nPDF] = myNormalDistribution1D(nptos)
% $$$ function [nPDF]=myGaussDistribution(nptos)
% $$$ function [nPDF]=myGaussDistribution(xmin,xmax,nptos)
%********************************************
%Use: [nPDF] = myNormalDistribution1D(nptos)
%----------------------------------------------------------------------
%sigmaXci = sqrt(0.01*xmax);
%GaussX = (1/(sigmaX*sqrt(2*pi))) * exp(-(xrng-xmean).^2/(2*sigmaX^2));
%----------------------------------------------------------------------
%********************************************

%%%%%%%%%%%%%%%%
%GAUSSIAN CURVE:
%%%%%%%%%%%%%%%%
%................
xmin=0;
xmax=10;
%................
dx=(xmax-xmin)/(nptos-1);
x=[xmin:dx:xmax];
xm=mean(x);
%..........................
sigma=sqrt(0.2*xmax); %sigma^2=1;
%..........................
beta=1/(sigma*sqrt(2*pi)); %ojo: hay que poner el parentesis al denominador!!
%..........................
PDF=beta*exp(-(((x-xm).^2)/(2*sigma^2)));
%..........................
PDF002=beta*exp(-0.5*((x-xm)/sigma).^2);
%..........................

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZED GAUSSIAN CURVE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%
nPDF = PDF/max(PDF);
nPDF002 = PDF002/max(PDF002);

%....................
% $$$ figure(100)
% $$$ subplot(2,2,1)
% $$$ plot(x,PDF,'r*')
% $$$ axis([-inf +inf, 0 1])
% $$$ subplot(2,2,2)
% $$$ plot(x,PDF002,'r*')
% $$$ axis([-inf +inf, 0 1])
% $$$ subplot(2,2,3)
% $$$ plot(x,nPDF,'r*')
% $$$ axis([-inf +inf, 0 1])
% $$$ subplot(2,2,4)
% $$$ plot(x,nPDF002,'r*')
% $$$ axis([-inf +inf, 0 1])
%....................
