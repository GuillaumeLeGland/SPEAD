function [YI]=myYearConversion365to360days(Y)

%................
ndays1=365;
ndays2=360;
%................
dx=ndays1/ndays2;
X=[1:ndays1];
XI=[dx:dx:ndays1];
%................
[YI]=interp1(X,Y,XI);
%................


