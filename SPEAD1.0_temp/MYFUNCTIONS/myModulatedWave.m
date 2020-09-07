function [Y]=myModulatedWave(dt,tmax,ndays)

%...........................
t=[dt:dt:tmax];
%...........................
% $$$ A0 = 1.0;
% $$$ C0 = 1.0;
% $$$ M0 = 0.5;
%...........................
A0 = 1.0;
C0 = 0.2;
M0 = 0.5;
%...........................
Tc=2*ndays;
Tm=1*ndays;
%...........................
fc = 1/10; %every 10 years modulation.
wc=(2*pi/Tc);
wm=(2*pi/Tm)*fc;
%...........................
% $$$ % $$$ phic=0*Tc;
% $$$ % $$$ phim=0*Tm;
% $$$ 
% $$$ phic=(3/4)*Tc;
% $$$ phim=(3/4)*Tm;
% $$$ 
% $$$ C = C0*(+sin(wc*t+phic));
% $$$ M = M0*(-cos(wm*t+phim));
%...........................
phic=0*Tc;
phim=0*Tm;

% $$$ phic=(3/4)*Tc;
% $$$ phim=(3/4)*Tm;

tc = t-phic;
tm = t-phim;

C = C0*(+sin(wc*tc));
M = M0*(-cos(wm*tm));
%...........................
Y = (A0 + M).*C;
%...........................
% $$$ figure(1)
% $$$ subplot(2,2,1)
% $$$ plot(C)
% $$$ grid on
% $$$ subplot(2,2,2)
% $$$ plot(M)
% $$$ grid on
% $$$ subplot(2,2,3)
% $$$ plot(Y)
% $$$ grid on
%...........................
