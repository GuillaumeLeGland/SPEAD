%%function [FX,FY,sigmax,sigmay,xaxis,yaxis,xm,ym,dx,dy] = jamstecrest_derivatives2D_gaussianfxy()

%%%%%%%%%%%%%%%%%%%%%%%%%%
%GAUSSIAN DISTRIBUTION 2D:
%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------------------
%<http://mathworld.wolfram.com/GaussianFunction.html>
%-----------------------------------------------------------------------------------
% <http://en.wikipedia.org/wiki/Gaussian_function>
%
% fxy = Amax * exp(-(a*(x-xm)^2 + c*(y-ym)^2 + 2*b*(x-xm)*(y-ym)))
%
% a = (cos(teta)^2 / 2*sigmax^2) + (sin(teta)^2 / 2*sigmay^2)
% c = (sin(teta)^2 / 2*sigmax^2) + (cos(teta)^2 / 2*sigmay^2)
% b = -(sin(2*teta) / (4*sigmax^2)) + (sin(2*teta) / (4*sigmay^2))
% teta = pi/3; 
%-----------------------------------------------------------------------------------
%===================================================================================
%...................................................................................
dx = 0.1; 
xmin = -10.0;
xmax = +10.0;
xaxis = [xmin:dx:xmax]; 
xm = mean(xaxis);
%...................................................................................
dy = 0.2; 
ymin = -10.0;
ymax = +10.0;
yaxis = [ymin:dy:ymax]; 
ym = mean(yaxis);
%...................................................................................
[XAXIS,YAXIS] = meshgrid(xaxis,yaxis);
%...................................................................................
msize = length(xaxis);
nsize = length(yaxis);
%...................................................................................
%===================================================================================
%...................................................................................
% $$$ xm = log(2.0);
% $$$ ym = ssti;
%...................................................................................
% $$$ sigmax = (xmax - xmin) / 3; 
% $$$ sigmay = (ymax - ymin) / 3; 
%...................................................................................
% $$$ sigmax = (xmax - xmin) / 4; 
% $$$ sigmay = (ymax - ymin) / 2; 
%...................................................................................
% $$$ sigmax = (xmax - xmin) / 10; 
% $$$ sigmay = (ymax - ymin) / 10; 
%...................................................................................
sigmax = (xmax - xmin) / 10; %USAR ESTAS!!!
sigmay = (ymax - ymin) / 20; 
%...................................................................................
FXY001 = 1/(2*pi*sigmax*sigmay) * exp(-( (XAXIS - xm).^2 / (2*sigmax^2) + (YAXIS - ym).^2 / (2*sigmay^2) ));
%...................................................................................
%===================================================================================
%...................................................................................
FX = 1/sqrt(2*pi*sigmax*sigmax) * exp(- (XAXIS - xm).^2 / (2*sigmax^2) );
FY = 1/sqrt(2*pi*sigmay*sigmay) * exp(- (YAXIS - ym).^2 / (2*sigmay^2) );
%...................................................................................
FXY002 = FX .* FY;
%...................................................................................
%===================================================================================
%...................................................................................
%%Amax = 1.0;
Amax = 1.0/(2*pi*sigmax*sigmay);
%...................................................................................
teta = pi/2; 
%%teta = pi/3; 
%...................................................................................
ax =   (cos(teta)^2 / 2*sigmax^2) + (sin(teta)^2 / 2*sigmay^2); %Original from wikipedia.
ay =   (sin(teta)^2 / 2*sigmax^2) + (cos(teta)^2 / 2*sigmay^2);
axy = -(sin(2*teta) / 4*sigmax^2) + (sin(2*teta) / 4*sigmay^2);
%...................................................................................
ax = ax/(sigmax*sigmay)^2; %I need this reescaling to make it work.
ay = ay/(sigmax*sigmay)^2; 
axy = axy/1;
%...................................................................................
FXY003 = Amax * exp(-(ax*(XAXIS-xm).^2 + ay*(YAXIS-ym).^2 + 2*axy*(XAXIS-xm).*(YAXIS-ym)));
%...................................................................................
%===================================================================================
%...................................................................................
% $$$ FXY = FXY001;
FXY = FXY002;
%...................................................................................
figure(60)
subplot(2,2,1)
imagesc(xaxis,yaxis,FXY001)
colorbar('vertical')
grid on
subplot(2,2,2)
imagesc(xaxis,yaxis,FXY002)
colorbar('vertical')
grid on
subplot(2,2,3)
imagesc(xaxis,yaxis,FXY003)
colorbar('vertical')
grid on
%...................................................................................
%%return
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GAUSSIAN-POWER DISTRIBUTION 2D:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
alfa = 1.5;
%...................................................................................
HX = (FX.^alfa);
HY = (FY.^alfa);
%...................................................................................
HXY001 = (FXY.^alfa); 
HXY002 = (FX.^alfa) .* (FY.^alfa);
%...................................................................................
HXY = HXY002; 
%...................................................................................
FXalfa = HX;
FYalfa = HY;
FXYalfa = HXY; 
%...................................................................................
%===================================================================================
%NUMERICAL DERIVATIVES OF H(x,y):
%...................................................................................
nth1 = 1; %First derivative (gradient)
nth2 = 2; %Second derivative (laplacian)
%...................................................................................
ndimY = 1; %diff between rows (ie. gradient in y-direction or columnwise)
ndimX = 2; %diff between cols (ie. gradient in x-direction or rowes-wise)
%...................................................................................
delHXYdx   = diff(HXY,nth1,ndimX)   /dx;
delHXYdy   = diff(HXY,nth1,ndimY)   /dy;
%...................................................................................
% $$$ delHXYdx = [delHXYdx,ones(nsize,1)*nan]; 
% $$$ delHXYdy = [delHXYdy;ones(1,msize)*nan]; 
%...................................................................................
del2HXYdxx = diff(delHXYdx,nth1,ndimX)   /dx;
del2HXYdyy = diff(delHXYdy,nth1,ndimY)   /dy;
%...................................................................................
% $$$ del2HXYdxxBis = diff(HXY,nth2,ndimX)   /dx; %SALE MAL!!! (by a factor dx)
% $$$ del2HXYdyyBis = diff(HXY,nth2,ndimY)   /dy; %SALE MAL!!! (by a factor dy)
%...................................................................................
del2HXYdxy = diff(delHXYdy,nth1,ndimX)/dx;
del2HXYdyx = diff(delHXYdx,nth1,ndimY)/dy;
%...................................................................................
[graHXYdx  ,graHXYdy]   = gradient(HXY     ,dx,dy);
[gra2HXYdxx,gra2HXYdyx] = gradient(graHXYdx,dx,dy);
[gra2HXYdxy,gra2HXYdyy] = gradient(graHXYdy,dx,dy);
%...................................................................................
%===================================================================================
%ANALYTICAL DERIVATIVES OF H(x,y):
%...................................................................................
sx = - (xaxis - xm) / sigmax^2; 
sy = - (yaxis - ym) / sigmay^2; 
%...................................................................................
[SX, junkarg] = meshgrid(sx,yaxis);
[junkarg, SY] = meshgrid(xaxis,sy);
%...................................................................................
dHXYdx = alfa * HXY .* SX; 
dHXYdy = alfa * HXY .* SY; 
%...................................................................................
d2HXYdxx = alfa * HXY .* ( alfa * (SX).^2 - (1/sigmax^2) );
d2HXYdyy = alfa * HXY .* ( alfa * (SY).^2 - (1/sigmay^2) );
%...................................................................................
d2HXYdxy001 = (alfa * FX.^alfa .* SX) .* (alfa * FY.^alfa .* SY); %OKAY.
d2HXYdxy002 = (alfa * HX .* SX) .* (alfa * HY .* SY);
d2HXYdxy003 = (alfa^2) * (HX .* HY) .* (SX .* SY);
d2HXYdxy004 = (alfa^2) * (HXY) .* (SX .* SY);
%...................................................................................
d2HXYdxy = d2HXYdxy002;
d2HXYdyx = d2HXYdxy; 
%...................................................................................
% $$$ figure(1000)
% $$$ % $$$ subplot(2,2,1)
% $$$ % $$$ plot(del2HXYdxx(:),d2HXYdxx(:),'*')
% $$$ % $$$ grid on 
% $$$ % $$$ subplot(2,2,2)
% $$$ % $$$ plot(del2HXYdyy(:),d2HXYdyy(:),'*')
% $$$ % $$$ grid on 
% $$$ subplot(2,2,3)
% $$$ plot(gra2HXYdxx(:),d2HXYdxx(:),'*')
% $$$ axis([-10d-3 +10d-3, -10d-3 +10d-3])
% $$$ grid on 
% $$$ subplot(2,2,4)
% $$$ plot(gra2HXYdyy(:),d2HXYdyy(:),'*')
% $$$ axis([-0.05 +0.05, -0.05 +0.05])
% $$$ grid on 
% $$$ %%return
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NUMERICAL INTEGRATION OF GASSIAN DISTRUBUTION 2D:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
ejex = 2;
ejey = 1;
%...................................................................................
sumFXdx = sum(FX*dx,ejex); %Should be equal to sqrt(1).
sumFYdy = sum(FY*dy,ejey); %Should be equal to sqrt(1).
%...................................................................................
sumFXalfadx = sum(FXalfa*dx,ejex); %Should be < than sqrt(1) and varies with sigmax.
sumFYalfady = sum(FYalfa*dy,ejey); %Should be < than sqrt(1) and varies with sigmay.
%...................................................................................
sumFXYdxdy = sum(sum(FXY*dx,ejex)*dy,ejey); %Should be equal to 1.
%...................................................................................
sumFXYalfadxdy = sum(sum(FXYalfa*dx,ejex)*dy,ejey); %Should be < than 1 and varies with sigmay.
%...................................................................................
%===================================================================================
%...................................................................................
Fx = FX(1,:);
Fy = FY(:,1);
%...................................................................................
Fx = FX(1,:);
Fy = FY(:,1);
%...................................................................................
Fxalfa = FXalfa(1,:);
Fyalfa = FYalfa(:,1);
%...................................................................................
%===================================================================================
%...................................................................................
sumFxdx = sum(Fx*dx); %Should be equal to sqrt(1).
sumFydy = sum(Fy*dy); %Should be equal to sqrt(1).
%...................................................................................
sumFxalfadx = sum(Fxalfa*dx); %Should be < than sqrt(1) and varies with sigmax.
sumFyalfady = sum(Fyalfa*dy); %Should be < than sqrt(1) and varies with sigmay.
%...................................................................................
%===================================================================================
%...................................................................................
%ANALYTICAL INTEGRATION OF PHYTOPLANKTON GASSIAN DISTRUBUTION:
%...................................................................................
sigmaxStar = sigmax / sqrt(alfa);
sigmayStar = sigmay / sqrt(alfa);
%...................................................................................
intFxdx = 1.0 
intFydy = 1.0 
%...................................................................................
intFxalfadx001 = (1.0 / (sigmax*sqrt(2*pi)))^alfa * (sigmaxStar*sqrt(2*pi)); %OKAY
intFyalfady001 = (1.0 / (sigmay*sqrt(2*pi)))^alfa * (sigmayStar*sqrt(2*pi));
%...................................................................................
intFxalfadx002 = (1.0 / sqrt(alfa)) * (sigmax*sqrt(2*pi))^(1-alfa); %OKAY
intFyalfady002 = (1.0 / sqrt(alfa)) * (sigmay*sqrt(2*pi))^(1-alfa);
%...................................................................................
intFxalfadx = intFxalfadx001
intFyalfady = intFyalfady001
%...................................................................................
intFxydxdy = intFxdx * intFxdx %OKAY
%...................................................................................
intFxyalfadxdy001 = intFxalfadx * intFxalfadx; %WRONG!!!!
intFxyalfadxdy002 = (1.0 / ((sigmax*sigmay)*(2*pi)))^alfa * (sigmaxStar*sigmayStar)*(2*pi); %OKAY
intFxyalfadxdy003 = (1/alfa) * ((sigmax*sigmay)*(2*pi))^(1-alfa); %OKAY
%...................................................................................
intFxyalfadxdy = intFxyalfadxdy003; 
%...................................................................................
Acte = intFxyalfadxdy; 
%...................................................................................
%===================================================================================
%RATIOS NUMERICAL / ANALYTICAL INTEGRATION:
%...................................................................................
ratioFxdx = sumFxdx / intFxdx
ratioFydy = sumFydy / intFydy
%...................................................................................
ratioFxalfadx = sumFxalfadx / intFxalfadx
ratioFyalfady = sumFyalfady / intFyalfady
%...................................................................................
ratioFxydxdy = sumFXYdxdy / intFxydxdy
ratioFxyalfadxdy = sumFXYalfadxdy / intFxyalfadxdy 
%...................................................................................
return

%%%%%%%
%PLOTS:
%%%%%%%
%===================================================================================
%...................................................................................
figure(20)
subplot(2,3,1)
imagesc(xaxis,yaxis,delHXYdx)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('dh(x,y) / dx')
grid on
subplot(2,3,2)
imagesc(xaxis,yaxis,del2HXYdxx)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dxx')
grid on
subplot(2,3,3)
imagesc(xaxis,yaxis,del2HXYdxy)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dxy')
grid on
subplot(2,3,4)
imagesc(xaxis,yaxis,delHXYdy)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('dh(x,y) / dy')
grid on
subplot(2,3,5)
imagesc(xaxis,yaxis,del2HXYdyy)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dyy')
grid on
subplot(2,3,6)
imagesc(xaxis,yaxis,del2HXYdyx)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dyx')
grid on
%...................................................................................
figure(30)
subplot(2,3,1)
imagesc(xaxis,yaxis,graHXYdx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('dh(x,y) / dx')
colorbar('horizontal')
grid on
subplot(2,3,2)
imagesc(xaxis,yaxis,gra2HXYdxx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dxx')
colorbar('horizontal')
grid on
subplot(2,3,3)
imagesc(xaxis,yaxis,gra2HXYdxy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dxy')
colorbar('horizontal')
grid on
subplot(2,3,4)
imagesc(xaxis,yaxis,graHXYdy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('dh(x,y) / dy')
colorbar('horizontal')
grid on
subplot(2,3,5)
imagesc(xaxis,yaxis,gra2HXYdyy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dyy')
colorbar('horizontal')
grid on
subplot(2,3,6)
imagesc(xaxis,yaxis,gra2HXYdyx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dyx')
colorbar('horizontal')
grid on
%...................................................................................
figure(40)
subplot(2,3,1)
imagesc(xaxis,yaxis,dHXYdx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('dh(x,y) / dx')
colorbar('horizontal')
grid on
subplot(2,3,2)
imagesc(xaxis,yaxis,d2HXYdxx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dxx')
colorbar('horizontal')
grid on
subplot(2,3,3)
imagesc(xaxis,yaxis,d2HXYdxy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dxy')
colorbar('horizontal')
grid on
subplot(2,3,4)
imagesc(xaxis,yaxis,dHXYdy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('dh(x,y) / dy')
colorbar('horizontal')
grid on
subplot(2,3,5)
imagesc(xaxis,yaxis,d2HXYdyy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dyy')
colorbar('horizontal')
grid on
subplot(2,3,6)
imagesc(xaxis,yaxis,d2HXYdyx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2h(x,y) / dyx')
colorbar('horizontal')
grid on
%...................................................................................
%%return
%===================================================================================

%***********************************************************************************
return
