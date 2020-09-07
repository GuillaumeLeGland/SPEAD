function [PDF] = myGaussPDF(Xrange,meanx,sigmax)

% $$$ Xrange = [-4:0.1:+4];
% $$$ meanx = 0;
% $$$ sigmax = 1;

%%%%%%%
%INPUT:
%%%%%%%
X = Xrange - nanmean(Xrange);
xm = meanx - nanmean(Xrange);
xsig = sigmax / std(Xrange);

w = [xm,xsig];

Y = normpdf(X,0,1);

plot(X,Y)
grid on

%%%%%%%%
%OUTPUT:
%%%%%%%%
PDF = Y;

