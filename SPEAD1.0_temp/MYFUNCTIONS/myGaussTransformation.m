function [Xgauss]=myGaussTransformation(Xdata)

%*************************************
%[Xgauss]=myGaussTransformation(Xdata)
%*************************************

%.............................
Xdata = Xdata(:);
msize = length(Xdata);
%.............................
%%npower = 0.5;
npower = 0.25;
%.............................
% $$$ for iday = 1:msize
% $$$     Xgauss(iday) = (Xdata(iday)^npower - 1.0) / npower*(geomean(Xdata)^(npower - 1.0));
% $$$ end
%.............................
%%gmeanXdata = geomean(Xdata);
%.............................
gmeanXdata = nangeomean(Xdata);
%.............................
Xgauss1 = (Xdata.^npower - 1.0) / npower*(gmeanXdata^(npower - 1.0));
%.............................
A0 = (0^npower - 1.0) / npower; %To force "zero" to be the minimum value.
%.............................
Xgauss2 = A0*(-1) + (Xdata.^npower - 1.0) / npower; %BoxCox.
%.............................

%%%%%%%
%PLOTS:
%%%%%%%
%.............................
figure(1)
subplot(2,2,1)
hist(Xdata(:),100)
subplot(2,2,3)
hist(Xgauss1(:),100)
subplot(2,2,4)
hist(Xgauss2(:),100)
%.............................

%%%%%%%%
%OUTPUT:
%%%%%%%%
%%Xgauss = Xgauss1;
Xgauss = Xgauss2;
