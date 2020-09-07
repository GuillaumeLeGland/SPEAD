close all
clear all
addpath ~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYFUNCTIONS/
GenPath=genpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYTOOLBOX/');
addpath(GenPath)

%ANALYTICAL DERIVATIVES OF GAUSSIAN FUNCTION:
%-----------------------------------------------------------------------------------
%<http://statistics.about.com/od/Mathstat/a/Inflection-Points-Of-The-Probability-Density-Function-Of-A-Normal-Distribution.htm>
%-----------------------------------------------------------------------------------
%===================================================================================
%...................................................................................
dx = 0.1; 
xmin = -10.0;
xmax = +10.0;
x = [xmin:dx:xmax]; 
xm = mean(x);
sigmax = 2.0;
%...................................................................................
a = -1/(2*sigmax^2); 
b = -(2*a)*xm; 
c = b^2 / (4*a) + log(-a/pi) / 2; 
%...................................................................................
fx000 = exp((a*x.^2) + (b*x) + c); %General expression.
fx001 = (1.0 / (sigmax * sqrt(2*pi))) * exp( -(1/2)*((x - xm)/sigmax).^2 ); %Okay.
fx002 = (1.0 / (sigmax * sqrt(2*pi))) * exp( -(x - xm).^2 / (2*sigmax^2) ); %Okay.
%...................................................................................
fx = fx002; 
%...................................................................................
dfxdx001 = -(x - xm) / (sigmax^3*sqrt(2*pi) ) .* exp(-(x - xm).^2 / (2*sigmax^2)); 
dfxdx002 = -(x - xm) .* (fx  / sigmax^2); 
dfxdx003 =  fx .* (-(x - xm) / sigmax^2); 
%...................................................................................
dfxdx  = dfxdx003; 
%...................................................................................
d2fxdx001 = - (fx./sigmax^2) - (x - xm) .* (dfxdx/sigmax^2); 
d2fxdx002 = - (fx./sigmax^2) + (x - xm).^2 .* (fx/(sigmax^4)); 
%...................................................................................
d2fxdx = d2fxdx002; 
%...................................................................................
figure(10)
subplot(2,3,1)
plot(x,fx)
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('size')
ylabel('PDF')
grid on
subplot(2,3,2)
plot(x,dfxdx001,'b-')
hold on
plot(x,dfxdx002,'r--')
hold off
grid on
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative first')
subplot(2,3,3)
plot(x,d2fxdx001,'b-')
hold on
plot(x,d2fxdx002,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative second')
grid on
%...................................................................................
%===================================================================================
%PHYTOPLANKTON BIOMASS GAUSSIAN DISTRIBUTION:
%-----------------------------------------------------------------------------------
%For alfa = 1.0 the g(x) becomes as constant (i.e. uniform distribution)
%For alfa = 2.0 the g(x) becomes as f(x) (i.e. original pure Gaussian)
%-----------------------------------------------------------------------------------
%...................................................................................
% $$$ Ctot = 1.0; 
% $$$ alfa = 1.0;
% $$$ Beta = 2.0;
%...................................................................................
% $$$ Ctot = 4.0; 
% $$$ alfa = 1.0;
% $$$ Beta = 2.0;
%...................................................................................
% $$$ Ctot = 1.0; 
% $$$ alfa = 2.0;
% $$$ Beta = 2.0;
%...................................................................................
% $$$ Ctot = 1.0; 
% $$$ alfa = 0.5;
% $$$ Beta = 2.0;
%...................................................................................
% $$$ Ctot = 3.0; 
% $$$ alfa = 1.5;
% $$$ Beta = 2.0;
%...................................................................................
Ctot = 3.0; 
alfa = 1.5;
Beta = 2.0;
%...................................................................................
% $$$ Ctot = 3.0; 
% $$$ alfa = 3.0;
% $$$ Beta = 2.0;
%...................................................................................
Px = Ctot * fx; 
%...................................................................................
Pxalfa = Px.^alfa; 
%...................................................................................
%===================================================================================
%...................................................................................
%NUMERICAL INTEGRATION OF PHYTOPLANKTON GASSIAN DISTRUBUTION:
sumPxdx = sum(Px*dx)
sumPxalfadx = sum(Pxalfa*dx)
%...................................................................................
%ANALYTICAL INTEGRATION OF PHYTOPLANKTON GASSIAN DISTRUBUTION:
intPxdx = Ctot 
intPxalfadx = (Ctot^alfa/sqrt(alfa)) * (sigmax*sqrt(2*pi))^(1-alfa) 
intPxalfadxBis = (Ctot/(sigmax*sqrt(2*pi)))^alfa * sigmax*(sqrt(2*pi)/sqrt(alfa)) 
%...................................................................................
%RATIOS NUMERICAL / ANALYTICAL INTEGRATION:
ratioPxdx = sumPxdx / intPxdx
ratioPxalfadx = sumPxalfadx / intPxalfadx
ratioPxalfadxBis = sumPxalfadx / intPxalfadxBis
%...................................................................................
%===================================================================================
%...................................................................................
%NUMERICAL FIRST DERIVATIVE OF PHYTOPLANKTON GAUSSIAN DISTRIBUTION:
delPxdx = diff(Px)/dx; 
delPxdx = [delPxdx,nan]; 
graPxdx = gradient(Px,dx); 
%...................................................................................
%NUMERICAL SECOND DERIVATIVE OF PHYTOPLANKTON GAUSSIAN DISTRIBUTION:
del2Pxdx = diff(delPxdx)/dx; 
del2Pxdx = [del2Pxdx,nan]; 
gra2Pxdx = gradient(graPxdx,dx); 
%...................................................................................
%ANALYTICAL FIRST DERIVATIVE OF PHYTOPLANKTON GAUSSIAN DISTRIBUTION:
dPxdx = -(x - xm) .* (Px / sigmax^2); 
%...................................................................................
%ANALYTICAL SECOND DERIVATIVE OF PHYTOPLANKTON GAUSSIAN DISTRIBUTION:
%...................................................................................
d2Pxdx001 = - (Px./sigmax^2) + (x - xm).^2 .* (Px/(sigmax^4)); 
%...................................................................................
d2Pxdx002 = Px .* ( ((x - xm).^2 / sigmax^4) - (1.0 / sigmax^2) ); 
%...................................................................................
d2Pxdx = d2Pxdx002;
%...................................................................................
%===================================================================================
%...................................................................................
%NUMERICAL FIRST DERIVATIVE OF PHYTOPLANKTON-ALFA GAUSSIAN DISTRIBUTION:
delPxalfadx = diff(Pxalfa)/dx; 
delPxalfadx = [delPxalfadx,nan]; 
graPxalfadx = gradient(Pxalfa,dx); 
%...................................................................................
%NUMERICAL SECOND DERIVATIVE OF PHYTOPLANKTON-ALFA GAUSSIAN DISTRIBUTION:
del2Pxalfadx = diff(delPxalfadx)/dx; 
del2Pxalfadx = [del2Pxalfadx,nan]; 
gra2Pxalfadx = gradient(graPxalfadx,dx); 
%...................................................................................
%ANALYTICAL FIRST DERIVATIVE OF PHYTOPLANKTON-ALFA GAUSSIAN DISTRIBUTION:
%................................................................................... 
dPxalfadx001 = Ctot^alfa * (sigmax*sqrt(2*pi))^(alfa-1) * sqrt(alfa) * (-(x - xm)./sigmax^2) .* Px ; %Does NOT work!!!!
%...................................................................................
dPxalfadx002 = Ctot^alfa * (sigmax*sqrt(2*pi))^(alfa-1) * (1/sqrt(alfa)) * (-(x - xm)./(sigmax/sqrt(alfa))^2) .* Px ; %Does NOT work!!!!
%...................................................................................
dPxalfadx003 = Ctot^alfa * (sigmax*sqrt(2*pi))^(alfa-1) * (1/sqrt(alfa)) * dPxdx; %Does NOT work!!!!
%...................................................................................
dPxalfadx005 = alfa * (Ctot^alfa * fx.^alfa) .* (- (x - xm) / sigmax^2); %Okay!!!
%...................................................................................
dPxalfadx006 = alfa * (Px.^alfa) .* (- (x - xm) / sigmax^2); %Okay!!!
%...................................................................................
dPxalfadx = dPxalfadx006; 
%...................................................................................
%ANALYTICAL SECOND DERIVATIVE OF PHYTOPLANKTON-ALFA GAUSSIAN DISTRIBUTION:
%...................................................................................
d2Pxalfadx001 = alfa*Pxalfa .* (alfa*((x - xm).^2 / sigmax^4) - (1.0/sigmax^2)); 
%...................................................................................
d2Pxalfadx = d2Pxalfadx001; 
%...................................................................................
figure(200)
subplot(2,2,1)
plot(delPxalfadx,dPxalfadx,'*')
grid on
subplot(2,2,2)
plot(del2Pxalfadx,d2Pxalfadx,'*')
grid on
%...................................................................................
%===================================================================================
%...................................................................................
%GRAZING FUNCTIONAL RESPONSE KTW:
Z = 1.0;
gzmax = 1.0;
%%kgz = 0.001; 
kgz = 0.5; 
%...................................................................................
Gx001 = (Pxalfa ./ sum(Pxalfa*dx)) * (sum(Px*dx) / (sum(Px*dx) + kgz)) * (gzmax*Z); %Numerical.
Gx002 = (Pxalfa ./ intPxalfadx) * (intPxdx / (intPxdx + kgz)) * (gzmax*Z); %Analytical.
%...................................................................................
Gx = Gx001; 
%...................................................................................
%NUMERICAL FIRST DERIVATIVE OF GRAZING KTW FUNCTION:
delGxdx = diff(Gx)/dx; 
delGxdx = [delGxdx,nan]; 
graGxdx = gradient(Gx,dx); 
%...................................................................................
%NUMERICAL SECOND DERIVATIVE OF GRAZING KTW FUNCTION:
del2Gxdx = diff(delGxdx)/dx; 
del2Gxdx = [del2Gxdx,nan]; 
gra2Gxdx = gradient(graGxdx,dx); 
%...................................................................................
%ANALYTICAL FIRST DERIVATIVE OF GRAZING KTW FUNCTION:
A = (sigmax*sqrt(2*pi))^(1-alfa) * (1/sqrt(alfa)); 
G = (Ctot^Beta / (kgz^Beta + Ctot^Beta)) * (gzmax*Z); 
B = (G/A); 
dGxdx = alfa * fx.^alfa .* (- (x - xm) / sigmax^2) * B;
%...................................................................................
%ANALYTICAL SECOND DERIVATIVE OF GRAZING KTW FUNCTION:
d2Gxdx = alfa * fx.^alfa .* (alfa * ((x - xm).^2 / sigmax^4) - (1/sigmax^2)) * B;
%...................................................................................
%===================================================================================
%...................................................................................
%BIOMASS SPECIFIC [d-1] GRAZING FUNCTIONAL RESPONSE KTW:
%...................................................................................
gx001 = Gx./Px; 
gx002 = fx.^(alfa-1) * (B/Ctot); 
%...................................................................................
gx = gx001;
%...................................................................................
figure(30)
subplot(2,2,1)
plot(x,Gx,'b-')
hold on
plot(x,gx,'r-')
hold off
grid on
subplot(2,2,2)
plot(Gx,gx)
grid on
%...................................................................................
%NUMERICAL DERIVATIVES OF BIOMASS SPECIFIC GRAZING FUNCTIONAL RESPONSE KTW:
%...................................................................................
delgxdx = diffxy(x,gx,[],1); %First derivative.
del2gxdx = diffxy(x,gx,[],2); %Second derivative. 
del3gxdx = diffxy(x,gx,[],3); %Third derivative. 
del4gxdx = diffxy(x,gx,[],4); %Fourth derivative. 
%...................................................................................
%ANALYTICAL DERIVATIVES OF BIOMASS SPECIFIC GRAZING FUNCTIONAL RESPONSE KTW:
%...................................................................................
dgxdx001 = (alfa - 1) * fx.^(alfa - 1) .* (-(x - xm) / sigmax^2 ) * (B/Ctot); 
dgxdx002 = (1 - alfa) * fx.^(alfa - 1) .* ( (x - xm) / sigmax^2 ) * (B/Ctot); 
dgxdx003 = (1 - alfa) * fx.^(alfa - 1) .* ( (x - xm) / sigmax^2 ) * (G/(Ctot*A)); 
dgxdx004 = gx * (1 - alfa) .* ( (x - xm) / sigmax^2 ); 
dgxdx005 = gx * (alfa - 1) .* (-(x - xm) / sigmax^2 ); 
%...................................................................................
dgxdx = dgxdx005;
%...................................................................................
d2gxdx001 = gx * (alfa - 1) .* ( (alfa - 1) * (-(x - xm) / sigmax^2).^2 - (1/sigmax^2));
d2gxdx002 = gx * (alfa - 1) .* ( (alfa - 1) * ( (x - xm) / sigmax^2).^2 - (1/sigmax^2));
%...................................................................................
d2gxdx = d2gxdx001;
%...................................................................................
%===================================================================================
%...................................................................................
%EXPONENTIAL FUNCTION: 
e0 = 0.10;
ae = 0.25;
ex = e0 * exp(ae*x);
%...................................................................................
%NUMERICAL DERIVATIVES OF EXPONENTIAL FUNCTION:
%...................................................................................
delexdx = diffxy(x,ex,[],1); %First derivative.
%...................................................................................
del2exdx = diffxy(x,ex,[],2); %Second derivative. 
%...................................................................................
%ANALYTICAL DERIVATIVES OF EXPONENTIAL FUNCTION:
%...................................................................................
dexdx =   ae * ex;
%...................................................................................
d2exdx = (ae^2) * ex; 
%...................................................................................
%===================================================================================
%PHYTOPLANKTON GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
%...................................................................................
DIN = 1.0; 
%...................................................................................
mu0 = 0.4;
ks0 = 0.1;
amu = 0.20;
aks = 0.80;
%...................................................................................
% $$$ mu0 = 1.0;
% $$$ ks0 = 0.1;
% $$$ amu = 1.0;
% $$$ aks = 2.0;
%...................................................................................
mux = mu0 * exp(amu*x);
ksx = ks0 * exp(aks*x);
%...................................................................................
lx = (ksx ./ (ksx + DIN));
qx = (DIN ./ (DIN + ksx));
ux = mux .* qx; 
%...................................................................................
% $$$ plot(ux)
%...................................................................................
%NUMERICAL DERIVATIVES OF GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
%...................................................................................
del1uxdx = diffxy(x,ux,[],1); %First derivative.
del2uxdx = diffxy(x,ux,[],2); %Second derivative. 
del3uxdx = diffxy(x,ux,[],3); %Third derivative.
del4uxdx = diffxy(x,ux,[],4); %Fourth derivative. 
%...................................................................................
d1uxdxNum = del1uxdx;
d2uxdxNum = del2uxdx;
d3uxdxNum = del3uxdx;
d4uxdxNum = del4uxdx;
%...................................................................................
delqxdx  = diffxy(x,qx,[],1);
del2qxdx = diffxy(x,qx,[],2);
del3qxdx = diffxy(x,qx,[],3);
del4qxdx = diffxy(x,qx,[],4);
%...................................................................................
dellxdx  = diffxy(x,lx,[],1);
del2lxdx = diffxy(x,lx,[],2);
del3lxdx = diffxy(x,lx,[],3);
del4lxdx = diffxy(x,lx,[],4);
%...................................................................................
%ANALYTICAL DERIVATIVES OF GROWTH UPTAKE RATE MICHAELIS MENTEN FUNCTION:
%...................................................................................
d1qxdx = - qx .* (aks * lx);
%%d1lxdx =   qx .* (aks * lx);
d1lxdx = - d1qxdx; 
%...................................................................................
d2qxdx = -aks * (d1lxdx .* qx + d1qxdx .* lx); %Okay 
d2lxdx = - d2qxdx; 
%...................................................................................
d3qxdx = -aks * (d2lxdx .* qx + d2qxdx .* lx + 2 * d1lxdx .* d1qxdx); 
d3lxdx = - d3qxdx; 
%...................................................................................
d1muxdx  = amu * mux;
d1ksxdx  = aks * ksx;
%...................................................................................
d2muxdx = (amu^2) * mux;
d2ksxdx = (aks^2) * ksx;
%...................................................................................
d1uxdx001 = (d1muxdx .* qx) + (mux .* d1qxdx); %Okay
d1uxdx002 = (d1muxdx .* qx) - (mux .* qx .* (aks * lx));
d1uxdx003 = (amu * mux) .* qx - (mux .* qx .* (aks * lx));
d1uxdx004 = qx .* ( (amu * mux) - (mux .* (aks * lx)) );
d1uxdx005 = (mux .* qx) .* ( amu - (aks * lx) );
d1uxdx006 = ux .* ( amu - (aks * lx) );
d1uxdx007 = ux .* ( amu -  aks * lx  ); %Okay
%...................................................................................
d1uxdx = d1uxdx007;
%...................................................................................
d2uxdx001 = d1uxdx .* (amu - aks * lx) - ux .* (aks .* d1lxdx); %Okay
%...................................................................................
d2uxdx002 = ux .* ( (amu - aks * lx).^2 - aks^2 * lx .* qx ); %Okay
%...................................................................................
d2uxdx = d2uxdx002;
%...................................................................................
% $$$ d3uxdx = d2uxdx .* (amu - aks .* lx) - 2*(aks * d1uxdx .* d1lxdx) - ux .* d2lxdx; %WRONG
%...................................................................................
% $$$ d4uxdx = d3uxdx .* (amu - aks .* lx) - aks * d2uxdx .* d1lxdx - ...
% $$$ 	 2*aks*(d2uxdx .* d1lxdx + d1uxdx .* d2lxdx) - ...
% $$$ 	 (d1uxdx .* d2lxdx + ux .* d3lxdx); %WRONG
%...................................................................................
%===================================================================================
%...................................................................................
d3uxdx = d2uxdx .* (amu - aks .* lx) - 2*(aks * d1uxdx .* d1lxdx) - (aks * ux .* d2lxdx); %Okay
%...................................................................................
d4uxdx001 = d3uxdx .* (amu - aks .* lx) ...
    -   aks * (d2uxdx .* d1lxdx) ...
    - 2*aks * (d2uxdx .* d1lxdx + d1uxdx .* d2lxdx) ...
    -   aks * (d1uxdx .* d2lxdx + ux .* d3lxdx); %Okay
%...................................................................................
d4uxdx002 = d3uxdx .* (amu - aks .* lx) ...
    - 3*aks * (d2uxdx .* d1lxdx) ...
    - 3*aks * (d1uxdx .* d2lxdx) ...
    -   aks * (  ux   .* d3lxdx); %Okay
%...................................................................................
d4uxdx003 = d3uxdx .* (amu - aks .* lx) ...
    - 3*aks * (d2uxdx .* d1lxdx + d1uxdx .* d2lxdx) ...
    -   aks * (  ux   .* d3lxdx); %Okay
%...................................................................................
d4uxdx = d4uxdx003;
%...................................................................................
d1uxdxSer = d1uxdx;
d2uxdxSer = d2uxdx;
d3uxdxSer = d3uxdx;
d4uxdxSer = d4uxdx;
%...................................................................................
%===================================================================================
%LAN SMITH DERIVATIVES:
%...................................................................................
d1uxdx = (amu*mux.*qx) - aks*(ksx./(mux.*DIN)).*(mux.*qx).^2;
%...................................................................................
d2uxdx = (amu - 2*aks*((ksx.*mux.*qx)./(mux.*DIN))) .* d1uxdx ...
    - aks*(aks - amu)*(ksx./(mux.*DIN)).*(mux.*qx).^2;
%...................................................................................
d3uxdx = (amu - 2*aks*((ksx.*mux.*qx)./(mux.*DIN))) .* d2uxdx ...
    - 2*aks*(ksx./(mux.*DIN)) .* d1uxdx.^2 ...
    - 4*aks*(aks - amu)    .* ((ksx.*mux.*qx)./(mux.*DIN)) .* d1uxdx ...
    -   aks*(aks - amu).^2 .* ( ksx./(mux.*DIN)) .* ux.^2;
%...................................................................................
d4uxdx = (amu - 2*aks*(mux.*qx.*ksx./(mux.*DIN))) .* d3uxdx ...
    - 6*(aks*(aks-amu))*(mux.*qx.*ksx./(mux.*DIN)) .* d2uxdx ...
    - 6*(aks*(ksx./(mux.*DIN))) .* (d1uxdx .* d2uxdx) ...
    - 6*(aks*(aks-amu))*(ksx./(mux.*DIN)) .* d1uxdx.^2 ...
    - 6*(aks*(aks-amu).^2)*(mux.*qx.*ksx./(mux.*DIN)) .* d1uxdx ...
    - (aks*(aks-amu).^3)*(ksx./(mux.*DIN)) .* (mux.*qx).^2;
%...................................................................................
d1uxdxLan = d1uxdx;
d2uxdxLan = d2uxdx;
d3uxdxLan = d3uxdx;
d4uxdxLan = d4uxdx;
%...................................................................................
%===================================================================================
%SYMBOLIC DERIVATIVES USING MATLAB:
%...................................................................................
[d1uxdxSymbolic,d2uxdxSymbolic,d3uxdxSymbolic,d4uxdxSymbolic] = jamstecrest_derivatives_symbolic;
%...................................................................................
d1uxdxSym = eval(d1uxdxSymbolic); %Convert from symbolic to numerical values at double precision.
d2uxdxSym = eval(d2uxdxSymbolic);
d3uxdxSym = eval(d3uxdxSymbolic);
d4uxdxSym = eval(d4uxdxSymbolic);
%...................................................................................
% $$$ d1uxdxSym = (2*DIN*exp(x/5))./(25*(DIN + exp((4*x)/5)/10)) - (4*DIN*exp(x/5).*exp((4*x)/5))./(125*(DIN + exp((4*x)/5)/10).^2);
% $$$ %...................................................................................
% $$$ d2uxdxSym = ( 2*DIN*exp(x/5))./(125*(DIN + exp((4*x)/5)/10)) - ...
% $$$ 	    (24*DIN*exp(x/5).*exp((4*x)/5))./( 625*(DIN + exp((4*x)/5)/10).^2) + ...
% $$$ 	    (16*DIN*exp(x/5).*exp((8*x)/5))./(3125*(DIN + exp((4*x)/5)/10).^3);
% $$$ %...................................................................................
% $$$ d3uxdxSym = (  2*DIN*exp(x/5))./(625*(DIN + exp((4*x)/5)/10)) - ...
% $$$ 	    (124*DIN*exp(x/5).*exp((4*x)/5))./(3125*(DIN + exp((4*x)/5)/10).^2) + ...
% $$$ 	    ( 48*DIN*exp(x/5).*exp((8*x)/5))./(3125*(DIN + exp((4*x)/5)/10).^3) - ...
% $$$ 	    ( 96*DIN*exp(x/5).*exp((4*x)/5).*exp((8*x)/5))./(78125*(DIN + exp((4*x)/ 5)/10).^4);
% $$$ %...................................................................................
% $$$ d4uxdxSym = (   2*DIN*exp(x/5))./(3125*(DIN + exp((4*x)/5)/10)) - ...
% $$$ 	    ( 624*DIN*exp(x/5).*exp(( 4*x)/5))./(  15625*(DIN + exp((4*x)/5)/10).^2) + ...
% $$$ 	    (2656*DIN*exp(x/5).*exp(( 8*x)/5))./(  78125*(DIN + exp((4*x)/5)/10).^3) + ...
% $$$ 	    ( 768*DIN*exp(x/5).*exp((16*x)/5))./(1953125*(DIN + exp((4*x)/5)/10).^5) - ...
% $$$ 	    (2688*DIN*exp(x/5).*exp(( 4*x)/5).*exp((8*x)/5))./(390625*(DIN + exp((4*x)/5)/10).^4);
%...................................................................................
%===================================================================================
%...................................................................................
[Absdistd1uxdxSer,Pcndistd1uxdxSer] = jamstecrest_derivatives1D_comparison(d1uxdxSym,d1uxdxSer);
[Absdistd2uxdxSer,Pcndistd2uxdxSer] = jamstecrest_derivatives1D_comparison(d2uxdxSym,d2uxdxSer);
[Absdistd3uxdxSer,Pcndistd3uxdxSer] = jamstecrest_derivatives1D_comparison(d3uxdxSym,d3uxdxSer);
[Absdistd4uxdxSer,Pcndistd4uxdxSer] = jamstecrest_derivatives1D_comparison(d4uxdxSym,d4uxdxSer);
%...................................................................................
[Absdistd1uxdxLan,Pcndistd1uxdxLan] = jamstecrest_derivatives1D_comparison(d1uxdxSym,d1uxdxLan);
[Absdistd2uxdxLan,Pcndistd2uxdxLan] = jamstecrest_derivatives1D_comparison(d2uxdxSym,d2uxdxLan);
[Absdistd3uxdxLan,Pcndistd3uxdxLan] = jamstecrest_derivatives1D_comparison(d3uxdxSym,d3uxdxLan);
[Absdistd4uxdxLan,Pcndistd4uxdxLan] = jamstecrest_derivatives1D_comparison(d4uxdxSym,d4uxdxLan);
%...................................................................................
% $$$ [Absdistd1uxdxNum,Pcndistd1uxdxNum] = jamstecrest_derivatives1D_comparison(d1uxdxSym,d1uxdxNum);
% $$$ [Absdistd2uxdxNum,Pcndistd2uxdxNum] = jamstecrest_derivatives1D_comparison(d2uxdxSym,d2uxdxNum);
% $$$ [Absdistd3uxdxNum,Pcndistd3uxdxNum] = jamstecrest_derivatives1D_comparison(d3uxdxSym,d3uxdxNum);
% $$$ [Absdistd4uxdxNum,Pcndistd4uxdxNum] = jamstecrest_derivatives1D_comparison(d4uxdxSym,d4uxdxNum);
%...................................................................................
[Absdistd1uxdxSerLan,Pcndistd1uxdxSerLan] = jamstecrest_derivatives1D_comparison(d1uxdxLan,d1uxdxSer);
[Absdistd2uxdxSerLan,Pcndistd2uxdxSerLan] = jamstecrest_derivatives1D_comparison(d2uxdxLan,d2uxdxSer);
[Absdistd3uxdxSerLan,Pcndistd3uxdxSerLan] = jamstecrest_derivatives1D_comparison(d3uxdxLan,d3uxdxSer);
[Absdistd4uxdxSerLan,Pcndistd4uxdxSerLan] = jamstecrest_derivatives1D_comparison(d4uxdxLan,d4uxdxSer);
%...................................................................................
close all
%===================================================================================
%...................................................................................
u001Ser = [d1uxdx(:),d1uxdxSer(:),Absdistd1uxdxSer,Pcndistd1uxdxSer];%,pause
u002Ser = [d2uxdx(:),d2uxdxSer(:),Absdistd2uxdxSer,Pcndistd2uxdxSer];%,pause
u003Ser = [d3uxdx(:),d3uxdxSer(:),Absdistd3uxdxSer,Pcndistd3uxdxSer];%,pause
u004Ser = [d4uxdx(:),d4uxdxSer(:),Absdistd4uxdxSer,Pcndistd4uxdxSer];%,pause 
%...................................................................................
u001Lan = [d1uxdx(:),d1uxdxLan(:),Absdistd1uxdxLan,Pcndistd1uxdxLan];%,pause
u002Lan = [d2uxdx(:),d2uxdxLan(:),Absdistd2uxdxLan,Pcndistd2uxdxLan];%,pause
u003Lan = [d3uxdx(:),d3uxdxLan(:),Absdistd3uxdxLan,Pcndistd3uxdxLan];%,pause
u004Lan = [d4uxdx(:),d4uxdxLan(:),Absdistd4uxdxLan,Pcndistd4uxdxLan];%,pause 
%...................................................................................
% $$$ u001Num = [d1uxdx(:),d1uxdxNum(:),Absdistd1uxdxNum,Pcndistd1uxdxNum];%,pause
% $$$ u002Num = [d2uxdx(:),d2uxdxNum(:),Absdistd2uxdxNum,Pcndistd2uxdxNum];%,pause
% $$$ u003Num = [d3uxdx(:),d3uxdxNum(:),Absdistd3uxdxNum,Pcndistd3uxdxNum];%,pause
% $$$ u004Num = [d4uxdx(:),d4uxdxNum(:),Absdistd4uxdxNum,Pcndistd4uxdxNum];%,pause 
%...................................................................................
u001SerLan = [d1uxdxSer(:),d1uxdxLan(:),Absdistd1uxdxSerLan,Pcndistd1uxdxSerLan];%,pause
u002SerLan = [d2uxdxSer(:),d2uxdxLan(:),Absdistd2uxdxSerLan,Pcndistd2uxdxSerLan];%,pause
u003SerLan = [d3uxdxSer(:),d3uxdxLan(:),Absdistd3uxdxSerLan,Pcndistd3uxdxSerLan];%,pause
u004SerLan = [d4uxdxSer(:),d4uxdxLan(:),Absdistd4uxdxSerLan,Pcndistd4uxdxSerLan];%,pause 
%...................................................................................
%===================================================================================
%...................................................................................
% $$$ d1uxdx = d1uxdxSer; 
% $$$ d2uxdx = d2uxdxSer; 
% $$$ d3uxdx = d3uxdxSer; 
% $$$ d4uxdx = d4uxdxSer; 
%...................................................................................
% $$$ d1uxdx = d1uxdxLan; 
% $$$ d2uxdx = d2uxdxLan; 
% $$$ d3uxdx = d3uxdxLan; 
% $$$ d4uxdx = d4uxdxLan; 
%...................................................................................
d1uxdx = d1uxdxSym; 
d2uxdx = d2uxdxSym; 
d3uxdx = d3uxdxSym; 
d4uxdx = d4uxdxSym; 
%...................................................................................
%===================================================================================
%NUMERICAL HIGHER ORDER DERIVATIVES: 
%...................................................................................
del3Pxdx = diffxy(x,Px,[],3); %Third derivative.
del4Pxdx = diffxy(x,Px,[],4); %Fourth derivative. 
%...................................................................................
del3Pxalfadx = diffxy(x,Pxalfa,[],3); %Third derivative.
del4Pxalfadx = diffxy(x,Pxalfa,[],4); %Fourth derivative. 
%...................................................................................
del3Gxdx = diffxy(x,Gx,[],3); %Third derivative.
del4Gxdx = diffxy(x,Gx,[],4); %Fourth derivative. 
%...................................................................................
gra3Pxdx = gradient(gradient(gradient(Px,dx),dx),dx); %Third derivative. 
gra4Pxdx = gradient(gradient(gradient(gradient(Px,dx),dx),dx),dx); %Fourth derivative. 
%...................................................................................
gra3Pxalfadx = gradient(gradient(gradient(Pxalfa,dx),dx),dx); %Third derivative. 
gra4Pxalfadx = gradient(gradient(gradient(gradient(Pxalfa,dx),dx),dx),dx); %Fourth derivative. 
%...................................................................................
gra3Gxdx = gradient(gradient(gradient(Gx,dx),dx),dx); %Third derivative. 
gra4Gxdx = gradient(gradient(gradient(gradient(Gx,dx),dx),dx),dx); %Fourth derivative. 
%...................................................................................
%===================================================================================
%...................................................................................
gragxdx = gradient(gx,dx); %First derivative. 
gra2gxdx = gradient(gradient(gx,dx),dx); %Second derivative. 
gra3gxdx = gradient(gradient(gradient(gx,dx),dx),dx); %Third derivative. 
gra4gxdx = gradient(gradient(gradient(gradient(gx,dx),dx),dx),dx); %Fourth derivative. 
%...................................................................................
gra1uxdx  = gradient(ux,dx); %First derivative. 
gra2uxdx = gradient(gradient(ux,dx),dx); %Second derivative. 
gra3uxdx = gradient(gradient(gradient(ux,dx),dx),dx); %Third derivative. 
gra4uxdx = gradient(gradient(gradient(gradient(ux,dx),dx),dx),dx); %Fourth derivative. 
%...................................................................................
%===================================================================================
%...................................................................................
%NEW TRAIT Y-AXIS FOR THE OPTIMAL TEMPERATURE FOR GROWTH:
%...................................................................................
sst0 = 15;
sstj = 25;
%...................................................................................
Q10 = 1.0;
%%Q10 = 2.0;
%...................................................................................
dy = 0.1; 
ymin = -10.0;
ymax = +40.0;
y = [ymin:dy:ymax]; 
ym = mean(y);
sigmay = 2.0;
%...................................................................................
fy = exp(-(y - ym).^2 / (2*sigmay^2));
q10 = Q10 ^((sstj - sst0)/10); 
qy = fy * q10;
%...................................................................................
%NUMERICAL DERIVATIVES: 
delqydy  = diffxy(y,qy,[],1);
del2qydy = diffxy(y,qy,[],2);
del3qydy = diffxy(y,qy,[],3);
del4qydy = diffxy(y,qy,[],4);
%...................................................................................
%ANALYTICAL DERIVATIVES:
d1qydy  = fy .* (  -(y - ym)/sigmay^2) * q10;
d2qydy = fy .* ( (-(y - ym)/sigmay^2).^2 - (1/sigmay^2) ) * q10;
%...................................................................................

%===================================================================================
%PLOTS:
%===================================================================================
% $$$ %...................................................................................
% $$$ delgxdxBis = delGxdx ./Px; %Uncorrect way of doing it!!!!
% $$$ del2gxdxBis = del2Gxdx./Px;
% $$$ del3gxdxBis = del3Gxdx./Px;
% $$$ del4gxdxBis = del4Gxdx./Px;
% $$$ %...................................................................................
% $$$ figure(31)
% $$$ subplot(2,2,1)
% $$$ plot(delgxdx,delgxdxBis)
% $$$ xlabel('d(Gx / Px)dx (correct)')
% $$$ ylabel('d(Gx)dx / Px (uncorrect)')
% $$$ set(gca,'Xlim',[-0.50 +0.50])
% $$$ set(gca,'Ylim',[-0.50 +0.50])
% $$$ title('First derivative')
% $$$ grid on
% $$$ subplot(2,2,2)
% $$$ plot(del2gxdx,del2gxdxBis)
% $$$ xlabel('d(Gx / Px)dx (correct)')
% $$$ ylabel('d(Gx)dx / Px (uncorrect)')
% $$$ set(gca,'Xlim',[-0.50 +0.50])
% $$$ set(gca,'Ylim',[-0.50 +0.50])
% $$$ title('Second derivative')
% $$$ grid on
% $$$ subplot(2,2,3)
% $$$ plot(del3gxdx,del3gxdxBis)
% $$$ xlabel('d(Gx / Px)dx (correct)')
% $$$ ylabel('d(Gx)dx / Px (uncorrect)')
% $$$ set(gca,'Xlim',[-0.50 +0.50])
% $$$ set(gca,'Ylim',[-0.50 +0.50])
% $$$ title('Third derivative')
% $$$ grid on
% $$$ subplot(2,2,4)
% $$$ plot(del4gxdx,del4gxdxBis)
% $$$ xlabel('d(Gx / Px)dx (correct)')
% $$$ ylabel('d(Gx)dx / Px (uncorrect)')
% $$$ set(gca,'Xlim',[-0.50 +0.50])
% $$$ set(gca,'Ylim',[-0.50 +0.50])
% $$$ title('Fourth derivative')
% $$$ grid on
% $$$ %...................................................................................
%===================================================================================
figure(20)
%...................................................................................
subplot(5,5,1)
plot(x,Px)
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('size')
ylabel('PDF - Gaussian')
grid on
subplot(5,5,2)
plot(x,delPxdx,'b-')
hold on
plot(x,graPxdx,'k-')
hold on
plot(x,dPxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 1st')
grid on
subplot(5,5,3)
plot(x,del2Pxdx,'b-')
hold on
plot(x,gra2Pxdx,'k-')
hold on
plot(x,d2Pxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 2nd')
grid on
subplot(5,5,4)
plot(x,del3Pxdx,'b-')
hold on
plot(x,gra3Pxdx,'k-')
hold on
%%plot(x,d3Pxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 3rd')
grid on
subplot(5,5,5)
plot(x,del4Pxdx,'b-')
hold on
plot(x,gra4Pxdx,'k-')
hold on
%%plot(x,d4Pxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 4th')
grid on
%...................................................................................
subplot(5,5,1+5)
plot(x,Pxalfa)
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('size')
ylabel('PDF - Qswitch')
grid on
subplot(5,5,2+5)
plot(x,delPxalfadx,'b-')
hold on
plot(x,graPxalfadx,'k-')
hold on
plot(x,dPxalfadx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 1st')
grid on
subplot(5,5,3+5)
plot(x,del2Pxalfadx,'b-')
hold on
plot(x,gra2Pxalfadx,'k-')
hold on
plot(x,d2Pxalfadx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 2nd')
grid on
subplot(5,5,4+5)
plot(x,del3Pxalfadx,'b-')
hold on
plot(x,gra3Pxalfadx,'k-')
hold on
%%plot(x,d3Pxalfadx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 3th')
grid on
subplot(5,5,5+5)
plot(x,del4Pxalfadx,'b-')
hold on
plot(x,gra4Pxalfadx,'k-')
hold on
%%plot(x,d4Pxalfadx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 4th')
grid on
%...................................................................................
subplot(5,5,1+10)
plot(x,Gx)
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('size')
ylabel('PDF - G(x)')
grid on
subplot(5,5,2+10)
plot(x,delGxdx,'b-')
hold on
plot(x,graGxdx,'k-')
hold on
plot(x,dGxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 1st')
grid on
subplot(5,5,3+10)
plot(x,del2Gxdx,'b-')
hold on
plot(x,gra2Gxdx,'k-')
hold on
plot(x,d2Gxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 2nd')
grid on
subplot(5,5,4+10)
plot(x,del3Gxdx,'b-')
hold on
plot(x,gra3Gxdx,'k-')
hold on
%%plot(x,d3Gxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 3th')
grid on
subplot(5,5,5+10)
plot(x,del4Gxdx,'b-')
hold on
plot(x,gra4Gxdx,'k-')
hold on
%%plot(x,d4Gxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 4th')
grid on
%...................................................................................
subplot(5,5,1+15)
plot(x,gx)
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('size')
ylabel('PDF - g(x)')
grid on
subplot(5,5,2+15)
plot(x,delgxdx,'b-')
hold on
plot(x,gragxdx,'k-')
hold on
plot(x,dgxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 1st')
grid on
subplot(5,5,3+15)
plot(x,del2gxdx,'b-')
hold on
plot(x,gra2gxdx,'k-')
hold on
plot(x,d2gxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 2nd')
grid on
subplot(5,5,4+15)
plot(x,del3gxdx,'b-')
hold on
plot(x,gra3gxdx,'k-')
hold on
%%plot(x,d3gxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 3th')
grid on
subplot(5,5,5+15)
plot(x,del4gxdx,'b-')
hold on
plot(x,gra4gxdx,'k-')
hold on
%%plot(x,d4gxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 4th')
grid on
%...................................................................................
subplot(5,5,1+20)
plot(x,ux)
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('size')
ylabel('PDF - u(x)')
grid on
subplot(5,5,2+20)
plot(x,del1uxdx,'b-')
hold on
plot(x,gra1uxdx,'k-')
hold on
plot(x,d1uxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 1st')
grid on
subplot(5,5,3+20)
plot(x,del2uxdx,'b-')
hold on
plot(x,gra2uxdx,'k-')
hold on
plot(x,d2uxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 2nd')
grid on
subplot(5,5,4+20)
plot(x,del3uxdx,'b-')
hold on
plot(x,gra3uxdx,'k-')
hold on
plot(x,d3uxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 3th')
grid on
subplot(5,5,5+20)
plot(x,del4uxdx,'b-')
hold on
plot(x,gra4uxdx,'k-')
hold on
plot(x,d4uxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('size')
ylabel('derivative 4th')
grid on
%...................................................................................
%===================================================================================
figure(25)
%...................................................................................
subplot(5,5,1)
plot(x,qx)
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('Temp Optima')
ylabel('PDF - q(x)')
grid on
subplot(5,5,2)
plot(x,delqxdx,'b-')
hold on
plot(x,d1qxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 1st')
grid on
subplot(5,5,3)
plot(x,del2qxdx,'b-')
hold on
plot(x,d2qxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 2nd')
grid on
subplot(5,5,4)
plot(x,del3qxdx,'b-')
hold on
plot(x,d3qxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 3th')
grid on
subplot(5,5,5)
plot(x,del4qxdx,'b-')
% $$$ hold on
% $$$ plot(x,d4qxdx,'r--')
% $$$ hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 4th')
grid on
%...................................................................................
subplot(5,5,1+05)
plot(x,lx)
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[0.00 1.00])
xlabel('Temp Optima')
ylabel('PDF - q(x)')
grid on
subplot(5,5,2+05)
plot(x,dellxdx,'b-')
hold on
plot(x,d1lxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 1st')
grid on
subplot(5,5,3+05)
plot(x,del2lxdx,'b-')
hold on
plot(x,d2lxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 2nd')
grid on
subplot(5,5,4+05)
plot(x,del3lxdx,'b-')
hold on
plot(x,d3lxdx,'r--')
hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 3th')
grid on
subplot(5,5,5+05)
plot(x,del4lxdx,'b-')
% $$$ hold on
% $$$ plot(x,d4lxdx,'r--')
% $$$ hold off
set(gca,'Xlim',[xmin xmax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 4th')
grid on
%...................................................................................
subplot(5,5,1+10)
plot(y,qy)
set(gca,'Xlim',[ymin ymax])
set(gca,'Ylim',[0.00 1.00])
xlabel('Temp Optima')
ylabel('PDF - q(y)')
grid on
subplot(5,5,2+10)
plot(y,delqydy,'b-')
hold on
% $$$ plot(y,graqydy,'k-')
% $$$ hold on
plot(y,d1qydy,'r--')
hold off
set(gca,'Xlim',[ymin ymax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 1st')
grid on
subplot(5,5,3+10)
plot(y,del2qydy,'b-')
hold on
% $$$ plot(y,gra2qydy,'k-')
% $$$ hold on
plot(y,d2qydy,'r--')
hold off
set(gca,'Xlim',[ymin ymax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 2nd')
grid on
subplot(5,5,4+10)
plot(y,del3qydy,'b-')
% $$$ hold on
% $$$ plot(y,gra3qydy,'k-')
% $$$ hold on
% $$$ %%plot(y,d3qydy,'r--')
% $$$ hold off
set(gca,'Xlim',[ymin ymax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 3th')
grid on
subplot(5,5,5+10)
plot(y,del4qydy,'b-')
% $$$ hold on
% $$$ plot(y,gra4qydy,'k-')
% $$$ hold on
% $$$ %%plot(y,d4qydy,'r--')
% $$$ hold off
set(gca,'Xlim',[ymin ymax])
set(gca,'Ylim',[-0.30 +0.30])
xlabel('Temp Optima')
ylabel('derivative 4th')
grid on
%...................................................................................
%===================================================================================
%***********************************************************************************
return

%%%%%%%%%%%%%%%%%%%%%%
%WEIBULL DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------
% <http://www.itl.nist.gov/div898/handbook/eda/section3/eda3668.htm> (Weibull distribution)
% <http://www.itl.nist.gov/div898/handbook/eda/section3/eda366b.htm> (Gamma distribution)
% <http://blog.minitab.com/blog/understanding-statistics/why-the-weibull-distribution-is-always-welcome> 
% <http://www.brighton-webs.co.uk/distributions/weibull2.aspx> 
%-----------------------------------------------------------------------------------
%...................................................................................
x = [xm:0.01:30.0]; 
%...................................................................................
delx = 10.0; %Location parameter.
rhox = 10.0; %Scale parameter (63.2 percentile of the data) 
%%ganmax = 1.0; %Shape parameter (if ganmax equal to 1.0 means that w(x) = exp(1/rhox); if ganmax equal to 3.5 means that w(x) = gaussian)
ganmax = 3.5; %Shape parameter (if ganmax equal to 1.0 means that w(x) = exp(1/rhox); if ganmax equal to 3.5 means that w(x) = gaussian)
%...................................................................................
wx001 = (ganmax/rhox) * (((x - delx)/rhox).^(ganmax - 1)) .* exp(-((x - delx)/rhox).^ganmax); 
%...................................................................................
figure(100)
plot(x,wx001) 
grid on
%...................................................................................
% <http://stackoverflow.com/questions/16965708/weibull-distribution-sample> 
Scale = 12.34;
Shape = 01.56;
[Mean,Variance] = wblstat(Scale,Shape)
figure(110)
ndata = 100;
mdata = 1;
subplot(1,2,1)
sample = wblrnd(Scale,Shape,ndata,mdata);    
histfit(sample,100,'wbl')
title('100 draws')
Scale = 12.34;
Shape = 01.56;
ndata = 1d5;
mdata = 1;
subplot(1,2,2)
sample = wblrnd(Scale,Shape,ndata,mdata);    
histfit(sample,100,'wbl')
title('100,000 draws')
%...................................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%
%LOG-NORMAL DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------
% <http://people.rit.edu/~jgvcqa/OtherFiles/WeibVsLN06_0523.pdf>
%-----------------------------------------------------------------------------------
% NOTE: 
% Sumatory of random quantities --> abs-Gaussian distribution 
% Products of random quantities --> log-Gaussian distribution 
%-----------------------------------------------------------------------------------
%NOTE: A random variable is lognormally distributed if the logarithm of the random variable is normally distributed.
%-----------------------------------------------------------------------------------
%...................................................................................
% $$$ dx = 0.1; 
% $$$ xmin = -10.0;
% $$$ xmax = +10.0;
% $$$ x = [xmin:dx:xmax] + abs(xmin); 
% $$$ xm = mean(x);
% $$$ sigmax = 2.0;
%...................................................................................
x = [0:0.01:20.0];
xm = 2.0;
sigmax = 2.0;
%...................................................................................
fx = (1.0 / (sigmax * sqrt(2*pi))) * exp( -(x - xm).^2 / (2*sigmax^2) ); %Okay.
%...................................................................................
logx = log(x);
logxm = log(xm); 
logsigmax = log(sigmax);
%...................................................................................
logfx001 = (1.0 ./ (x*logsigmax*sqrt(2*pi))) .* exp(- (logx - logxm).^2 / (2*logsigmax^2));
logfx002 = (1.0  / (logsigmax * sqrt(2*pi)))  * exp(- (logx - logxm).^2 / (2*logsigmax^2)); %Okay.
logfx003 = exp(fx); 
%...................................................................................
figure(200)
subplot(3,2,1)
plot(x,logfx001)
subplot(3,2,2)
plot(logx,logfx001)
subplot(3,2,3)
plot(x,logfx002)
subplot(3,2,4)
plot(logx,logfx002)
subplot(3,2,5)
plot(x,logfx003)
subplot(3,2,6)
plot(logx,logfx003)
%...................................................................................
%===================================================================================
%...................................................................................
% $$$ x = linspace(-1,2,20);
% $$$ y = exp(x);
%...................................................................................
x = [xmin:dx:xmax]; 
fx = fx; 
%...................................................................................
% $$$ dfx  = diffxy(x,fx);
% $$$ dfx2 = diffxy(x,dfx);  % Or, could use >> dy2 = diffxy(x,y,[],2);
% $$$ dfx3 = diffxy(x,dfx2);  % Or, could use >> dy2 = diffxy(x,y,[],2);
% $$$ dfx4 = diffxy(x,dfx3);  % Or, could use >> dy2 = diffxy(x,y,[],2);
%...................................................................................
dfx  = diffxy(x,fx);
dfx2 = diffxy(x,fx,[],2);  % Or, could use >> dy2 = diffxy(x,y,[],2);
dfx3 = diffxy(x,fx,[],3);  % Or, could use >> dy2 = diffxy(x,y,[],2);
dfx4 = diffxy(x,fx,[],4);  % Or, could use >> dy2 = diffxy(x,y,[],2);
%...................................................................................
Dfx  = gradient(fx)./gradient(x);
Dfx2 = gradient(Dfx)./gradient(x);
Dfx3 = gradient(Dfx2)./gradient(x);
Dfx4 = gradient(Dfx3)./gradient(x);
%...................................................................................
figure(300)
subplot(2,2,1)
plot(x,Dfx ,'-k',x,dfx ,'-r')
subplot(2,2,2)
plot(x,Dfx2,'-k',x,dfx2,'-r')
subplot(2,2,3)
plot(x,Dfx3,'-k',x,dfx3,'-r')
subplot(2,2,4)
plot(x,Dfx4,'-k',x,dfx4,'-r')
%...................................................................................
figure(310)
plot(x,(fx-dfx)./fx,'b*',x,(fx-dfx2)./fx,'b^')
hold on
plot(x,(fx-Dfx)./fx,'r*',x,(fx-Dfx2)./fx,'r^')
hold off
title('Relative error in derivative approximation')
legend('diffxy: dfx/dx','diffxy: d^2fx/dx^2','gradient: dfx/dx','gradient: d^2fx/dx^2')
hold off
%...................................................................................
%===================================================================================
