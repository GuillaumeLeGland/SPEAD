function [ysol]=mySolitonCurve(T,dt,tau,tmax,ymin,ymax)

%=================
%SINUSOIDAL CURVE:
%=================
%----------------------
%http://en.wikipedia.org/wiki/Sine_wave:
%
%x = A*sin(w*t + phi);
%----------------------
%...........
% $$$ T=365;
% $$$ tau=T/2;
% $$$ tmax=T*4;
% $$$ ymin=0; %[C]
% $$$ ymax=1;
%...........

t=[dt:dt:tmax];
tt=t-tau;
wrad=(2*pi/T);
A0=ymin;
A=(ymax-ymin)/2; %amplitud.
%..............
% $$$ fexp = sin(tt);
% $$$ ysol = A0 + A*(sin(wrad*tt)).* fexp;
%..............
% $$$ fexp = sin(tt);
% $$$ %%fexp = sin(tt)+1;
% $$$ ysol = A0 + A*( sin(wrad*tt) + 1) .* fexp;
%..............
% $$$ %%fexp = sin(wrad*tt.^2);
% $$$ fexp = sin((wrad*tt).^2);
% $$$ ysol = A0 + A*( sin(wrad*tt)) .* fexp;
%..............
% $$$ fexp = sin(wrad*tt).^2;
% $$$ ysol = A0 + A*( sin(wrad*tt)).*fexp + ymax;
%..............
% $$$ fexp = exp(sin(wrad*tt.^2));
% $$$ ysol = A0 + A*(sin(wrad*tt)).* fexp;
% $$$ %%ysol = A0 + A*( sin(wrad*tt) + 1) .* fexp;
%..............
% $$$ fexp = exp(sin(wrad*tt)); %USAR ESTA!!
% $$$ ysol = A0 + A*(sin(wrad*tt)).* fexp;
% $$$ %%ysol = A0 + A*( sin(wrad*tt) + 1) .* fexp;
%..............
% $$$ fgauss = exp(-(x-x0)/sigma^2)
%..............
% $$$ x=sin(wrad*tt);
% $$$ x0=0;
% $$$ sigma=1;
% $$$ fgauss = exp(-(x-x0).^2/sigma^2);
% $$$ ysol = fgauss;
%..............
% $$$ c=0;
% $$$ nptos=length(t);
% $$$ for i=1:nptos
% $$$     tti=tt(i);
% $$$     xi=wrad*tti;
% $$$     if c==0
% $$$ 	if floor(mod(xi,2*pi))==0
% $$$ 	    x0=xi;
% $$$ 	    c=c+1;
% $$$ 	end
% $$$     end
% $$$ % $$$     xi,x0
% $$$     sigma=1;
% $$$     fgauss(i) = exp(-(xi-x0).^2/sigma^2);
% $$$ end
% $$$ ysol = fgauss;
%..............
alfa=0.2;
beta=0.8;
nptos=length(t);
for i=1:nptos
    xi=tt(i);
    xa=floor(xi/pi)*pi + (pi/2);
% $$$     xb=xa+(2*pi);
    fgaussA = exp(-(xi-xa)^2/alfa);
% $$$     fgaussB = exp(-(xi-xb)^2/beta);
    fgaussB = 0;
    fgauss(i) =  fgaussA - fgaussB;
end
%%ysol = fgauss;
%%ysol = fgauss - ymax;
ysol = fgauss - mean(fgauss(:));  %USAR ESTA!!!!
%..............
Area=sum(ysol*dt)
%..............

figure(1)
plot(t,ysol,'.b-')

