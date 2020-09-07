%..............
dt=1;
ndays=360;
T=ndays;
tau=ndays/4;
tmax=T;
%..............
% $$$ ymin=0.1; %[m]
% $$$ ymax=1;
%..............
ymin=10; %[m]
ymax=150;
%..............
t=[dt:dt:tmax];
tt=t+tau;
wrad=(2*pi/T);
A0=ymin;
A=(ymax-ymin)/2; %amplitud.
%..............
sinwt = sin(wrad*tt);
%..............
ysin = A0 + A*( sinwt + 1);
%..............
%%xsin = (ysin/ymin).^2;
%%xsin = (ysin/ymin).^2 / ymax.^2;
xsin = (ysin/ymin).^2 / (ymax./ymin).^2; %USAR ESTA!!!!
%..............
npower=1.5;
%..............
ysin2 = ysin.*xsin;
ysin3 = A0 + 2*A*xsin;
ysin4 = A0 + 2*A*xsin.^2;
ysin5 = A0 + 2*A*xsin.^npower; %USAR ESTA!!!!
%..............

figure(100)
subplot(3,3,1)
plot(ysin*(-1))
grid on
subplot(3,3,2)
plot(ysin2*(-1))
grid on
subplot(3,3,3)
plot(ysin3*(-1))
grid on
subplot(3,3,4)
plot(ysin4*(-1))
grid on
subplot(3,3,5)
plot(ysin5*(-1))
grid on
subplot(3,3,9)
plot(xsin)
grid on
