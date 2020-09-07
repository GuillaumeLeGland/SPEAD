%%function [] = jamstecrest_derivatives2D_grazingKTW(FX,FY,sigmax,sigmay,xaxis,yaxis,xm,ym,dx,dy)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%GAUSSIAN DISTRIBUTION 2D:
%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
xmin = log(0.02); %USAR ESTOS!!!!
xmax = log(200.0);
dx = (xmax - xmin)/(256-1); 
xaxis = [xmin:dx:xmax]; 
%...................................................................................
ymin = 0; %log([degrees C])
ymax = 50; %log([degrees C])
dy = (ymax - ymin)/(128-1); 
yaxis  = [ymin:dy:ymax]; %Internal sst phy peak (maximum growth) [degress C]
%...................................................................................
xm = mean(xaxis);
ym = mean(yaxis);
%...................................................................................
sigmax = (xmax - xmin) / 10; %USAR ESTAS!!!
sigmay = (ymax - ymin) / 20; 
%...................................................................................
[XAXIS,YAXIS] = meshgrid(xaxis,yaxis);
%...................................................................................
FX = 1/sqrt(2*pi*sigmax*sigmax) * exp(- (XAXIS - xm).^2 / (2*sigmax^2) );
FY = 1/sqrt(2*pi*sigmay*sigmay) * exp(- (YAXIS - ym).^2 / (2*sigmay^2) );
%...................................................................................
FXY = FX .* FY;
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BIOMASS DISTRIBUTION PHYTOPLANKTON:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
%%Ctot = 1.0; 
Ctot = 9.0; 
alfa = 1.5;
Beta = 2.0;
%...................................................................................
msize = length(xaxis);
nsize = length(yaxis);
%...................................................................................
%===================================================================================
%...................................................................................
FXY = FX .* FY;
%...................................................................................
PXY001 = Ctot * FXY;
%...................................................................................
PXYalfa001 = (Ctot * FXY).^alfa;
%...................................................................................
%===================================================================================
%...................................................................................
Ctotx = sqrt(Ctot);
Ctoty = sqrt(Ctot);
%...................................................................................
PX = Ctotx * FX; 
PY = Ctoty * FY; 
%...................................................................................
PXalfa = PX.^alfa; 
PYalfa = PY.^alfa; 
%...................................................................................
PXY002 = PX .* PY;
%...................................................................................
PXYalfa002 = PXalfa .* PYalfa;
%...................................................................................
%===================================================================================
%...................................................................................
%%PXY = PXY001;
PXY = PXY002;
%...................................................................................
%%PXYalfa = PXYalfa001;
PXYalfa = PXYalfa002;
%...................................................................................
%===================================================================================
%...................................................................................
figure(70)
subplot(2,2,1)
imagesc(xaxis,yaxis,PXY001)
colorbar('vertical')
grid on
subplot(2,2,2)
imagesc(xaxis,yaxis,PXY002)
colorbar('vertical')
grid on
subplot(2,2,3)
imagesc(xaxis,yaxis,PXYalfa001)
colorbar('vertical')
grid on
subplot(2,2,4)
imagesc(xaxis,yaxis,PXYalfa002)
colorbar('vertical')
grid on
%...................................................................................
%%return
%===================================================================================
%NUMERICAL INTEGRATION OF PHYTOPLANKTON GASSIAN DISTRUBUTION 2D:
%...................................................................................
ejex = 2;
ejey = 1;
%...................................................................................
sumPXdx = sum(PX*dx,ejex); %Should be equal to sqrt(Ctot).
sumPYdy = sum(PY*dy,ejey); %Should be equal to sqrt(Ctot).
%...................................................................................
sumPXalfadx = sum(PXalfa*dx,ejex); %Should be < than sqrt(Ctot) and varies with sigmax.
sumPYalfady = sum(PYalfa*dy,ejey); %Should be < than sqrt(Ctot) and varies with sigmay.
%...................................................................................
sumPXYdxdy = sum(sum(PXY*dx,ejex)*dy,ejey); %Should be equal to Ctot.
sumPXYalfadxdy = sum(sum(PXYalfa*dx,ejex)*dy,ejey); %Should be < than Ctot and varies with sigmay.
%...................................................................................
%===================================================================================
%...................................................................................
Fx = FX(1,:);
Fy = FY(:,1);
%...................................................................................
Px = PX(1,:);
Py = PY(:,1);
%...................................................................................
Pxalfa = PXalfa(1,:);
Pyalfa = PYalfa(:,1);
%...................................................................................
%===================================================================================
%...................................................................................
sumPxdx = sum(Px*dx); %Should be equal to sqrt(Ctot).
sumPydy = sum(Py*dy); %Should be equal to sqrt(Ctot).
%...................................................................................
sumPxalfadx = sum(Pxalfa*dx); %Should be < than sqrt(Ctot) and varies with sigmax.
sumPyalfady = sum(Pyalfa*dy); %Should be < than sqrt(Ctot) and varies with sigmay.
%...................................................................................
% $$$ sumPxydxdy = sum(sum(Pxy*dx,ejex)*dy,ejey); %Should be equal to Ctot.
% $$$ sumPxyalfadxdy = sum(sum(Pxyalfa*dx,ejex)*dy,ejey); %Should be < than Ctot and varies with sigmay.
%...................................................................................
%===================================================================================
%...................................................................................
%ANALYTICAL INTEGRATION OF PHYTOPLANKTON GASSIAN DISTRUBUTION:
%...................................................................................
sigmaxStar = sigmax / sqrt(alfa);
sigmayStar = sigmay / sqrt(alfa);
%...................................................................................
intPxdx = Ctotx 
intPydy = Ctoty 
%...................................................................................
intPxalfadx001 = (Ctotx / (sigmax*sqrt(2*pi)))^alfa * (sigmaxStar*sqrt(2*pi)); %OKAY
intPyalfady001 = (Ctoty / (sigmay*sqrt(2*pi)))^alfa * (sigmayStar*sqrt(2*pi));
%...................................................................................
intPxalfadx002 = (Ctotx^alfa/sqrt(alfa)) * (sigmax*sqrt(2*pi))^(1-alfa); %OKAY
intPyalfady002 = (Ctoty^alfa/sqrt(alfa)) * (sigmay*sqrt(2*pi))^(1-alfa);
%...................................................................................
intPxalfadx = intPxalfadx001
intPyalfady = intPyalfady001
%...................................................................................
intPxydxdy = intPxdx * intPxdx 
%...................................................................................
Acte = (1/alfa) * ((sigmax*sigmay)*(2*pi))^(1-alfa);
%...................................................................................
intPxyalfadxdy001 = intPxalfadx * intPxalfadx; %WRONG!!!!
intPxyalfadxdy002 = (Ctot / ((sigmax*sigmay)*(2*pi)))^alfa * (sigmaxStar*sigmayStar)*(2*pi); %OKAY
intPxyalfadxdy003 = (Ctot^alfa) * (1/alfa) * ((sigmax*sigmay)*(2*pi))^(1-alfa); %OKAY
intPxyalfadxdy004 = (Ctot^alfa) * Acte; %OKAY
%...................................................................................
intPxyalfadxdy = intPxyalfadxdy004; 
%...................................................................................
%===================================================================================
%RATIOS NUMERICAL / ANALYTICAL INTEGRATION:
%...................................................................................
ratioPxdx = sumPxdx / intPxdx
ratioPydy = sumPydy / intPydy
%...................................................................................
ratioPxalfadx = sumPxalfadx / intPxalfadx
ratioPyalfady = sumPyalfady / intPyalfady
%...................................................................................
ratioPxydxdy = sumPXYdxdy / intPxydxdy
ratioPxyalfadxdy = sumPXYalfadxdy / intPxyalfadxdy 
%...................................................................................
return
%===================================================================================
%NUMERICAL DERIVATIVES OF PHYTOPLANKTON GAUSSIAN DISTRIBUTION 2D:
%...................................................................................
nth1 = 1; %First derivative (gradient)
nth2 = 2; %Second derivative (laplacian)
%...................................................................................
ndimY = 1; %diff between rows (ie. gradient in y-direction or columnwise)
ndimX = 2; %diff between cols (ie. gradient in x-direction or rowes-wise)
%...................................................................................
delPXYdx = diff(PXY,nth1,ndimX)/dx;
delPXYdy = diff(PXY,nth1,ndimY)/dy;
%...................................................................................
delPXYdx = [delPXYdx,ones(nsize,1)*nan]; 
delPXYdy = [delPXYdy;ones(1,msize)*nan]; 
%...................................................................................
del2PXYdxx = diff(delPXYdx,nth1,ndimX)/dx;
del2PXYdyy = diff(delPXYdy,nth1,ndimY)/dy;
%...................................................................................
del2PXYdxy = diff(delPXYdy,nth1,ndimX)/dx;
del2PXYdyx = diff(delPXYdx,nth1,ndimY)/dy;
%...................................................................................
[graPXYdx  ,graPXYdy]   = gradient(PXY     ,dx,dy);
[gra2PXYdxx,gra2PXYdyx] = gradient(graPXYdx,dx,dy);
[gra2PXYdxy,gra2PXYdyy] = gradient(graPXYdy,dx,dy);
%...................................................................................
%===================================================================================
%NUMERICAL DERIVATIVES OF PHYTOPLANKTON GAUSSIAN-POWER DISTRIBUTION 2D:
%...................................................................................
delPXYalfadx = diff(PXYalfa,nth1,ndimX)/dx;
delPXYalfady = diff(PXYalfa,nth1,ndimY)/dy;
%...................................................................................
delPXYalfadx = [delPXYalfadx,ones(nsize,1)*nan]; 
delPXYalfady = [delPXYalfady;ones(1,msize)*nan]; 
%...................................................................................
del2PXYalfadxx = diff(delPXYalfadx,nth1,ndimX)/dx;
del2PXYalfadyy = diff(delPXYalfady,nth1,ndimY)/dy;
%...................................................................................
del2PXYalfadxy = diff(delPXYalfady,nth1,ndimX)/dx;
del2PXYalfadyx = diff(delPXYalfadx,nth1,ndimY)/dy;
%...................................................................................
[graPXYalfadx  ,graPXYalfady]   = gradient(PXYalfa     ,dx,dy);
[gra2PXYalfadxx,gra2PXYalfadyx] = gradient(graPXYalfadx,dx,dy);
[gra2PXYalfadxy,gra2PXYalfadyy] = gradient(graPXYalfady,dx,dy);
%...................................................................................
%===================================================================================
%ANALYTICAL DERIVATIVES OF PHYTOPLANKTON GAUSSIAN DISTRIBUTION:
%...................................................................................
sx = - (xaxis - xm) / sigmax^2; 
sy = - (yaxis - ym) / sigmay^2; 
%...................................................................................
[SX, junkarg] = meshgrid(sx,yaxis);
[junkarg, SY] = meshgrid(xaxis,sy);
%...................................................................................
dPXYdx =  PXY .* SX; 
dPXYdy =  PXY .* SY; 
%...................................................................................
d2PXYdxx =  PXY .* ( (SX.^2) - (1/sigmax^2) );
d2PXYdyy =  PXY .* ( (SY.^2) - (1/sigmay^2) );
%...................................................................................
d2PXYdxy = (PXY) .* (SX .* SY);
d2PXYdyx = (PXY) .* (SX .* SY);
%...................................................................................
%===================================================================================
%ANALYTICAL DERIVATIVES OF PHYTOPLANKTON GAUSSIAN-POWER DISTRIBUTION:
%...................................................................................
dPXYalfadx = alfa * PXYalfa .* SX; 
dPXYalfady = alfa * PXYalfa .* SY; 
%...................................................................................
d2PXYalfadxx = alfa * PXYalfa .* ( alfa * (SX.^2) - (1/sigmax^2) );
d2PXYalfadyy = alfa * PXYalfa .* ( alfa * (SY.^2) - (1/sigmay^2) );
%...................................................................................
d2PXYalfadxy = (alfa^2) * (PXYalfa) .* (SX .* SY);
d2PXYalfadyx = (alfa^2) * (PXYalfa) .* (SX .* SY);
%...................................................................................
%===================================================================================
%ANALYTICAL DERIVATIVES OF GAUSSIAN DISTRIBUTION:
%...................................................................................
dFXYdx = FXY .* SX; 
dFXYdy = FXY .* SY; 
%...................................................................................
d2FXYdxx = FXY .* ( (SX.^2) - (1/sigmax^2) );
d2FXYdyy = FXY .* ( (SY.^2) - (1/sigmay^2) );
%...................................................................................
d2FXYdxy = (FXY) .* (SX .* SY);
d2FXYdyx = (FXY) .* (SX .* SY);
%...................................................................................
%===================================================================================
%ANALYTICAL DERIVATIVES OF GAUSSIAN-POWER DISTRIBUTION:
%...................................................................................
FXalfa = FX.^alfa;
FYalfa = FY.^alfa;
%...................................................................................
FXYalfa = FXY.^alfa;
%...................................................................................
dFXYalfadx = alfa * FXYalfa .* SX; 
dFXYalfady = alfa * FXYalfa .* SY; 
%...................................................................................
d2FXYalfadxx = alfa * FXYalfa .* ( alfa * (SX.^2) - (1/sigmax^2) );
d2FXYalfadyy = alfa * FXYalfa .* ( alfa * (SY.^2) - (1/sigmay^2) );
%...................................................................................
d2FXYalfadxy = (alfa^2) * (FXYalfa) .* (SX .* SY);
d2FXYalfadyx = (alfa^2) * (FXYalfa) .* (SX .* SY);
%...................................................................................
%===================================================================================
%...................................................................................
%...................................................................................
%...................................................................................
%===================================================================================

%%%%%%%
%PLOTS:
%%%%%%%
%===================================================================================
% $$$ %...................................................................................
% $$$ figure(20)
% $$$ subplot(2,3,1)
% $$$ imagesc(xaxis,yaxis,delPXYdx)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y) / dx')
% $$$ grid on
% $$$ subplot(2,3,2)
% $$$ imagesc(xaxis,yaxis,del2PXYdxx)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dxx')
% $$$ grid on
% $$$ subplot(2,3,3)
% $$$ imagesc(xaxis,yaxis,del2PXYdxy)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dxy')
% $$$ grid on
% $$$ subplot(2,3,4)
% $$$ imagesc(xaxis,yaxis,delPXYdy)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y) / dy')
% $$$ grid on
% $$$ subplot(2,3,5)
% $$$ imagesc(xaxis,yaxis,del2PXYdyy)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dyy')
% $$$ grid on
% $$$ subplot(2,3,6)
% $$$ imagesc(xaxis,yaxis,del2PXYdyx)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dyx')
% $$$ grid on
% $$$ %...................................................................................
% $$$ figure(30)
% $$$ subplot(2,3,1)
% $$$ imagesc(xaxis,yaxis,graPXYdx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y) / dx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,2)
% $$$ imagesc(xaxis,yaxis,gra2PXYdxx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dxx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,3)
% $$$ imagesc(xaxis,yaxis,gra2PXYdxy)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dxy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,4)
% $$$ imagesc(xaxis,yaxis,graPXYdy)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y) / dy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,5)
% $$$ imagesc(xaxis,yaxis,gra2PXYdyy)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dyy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,6)
% $$$ imagesc(xaxis,yaxis,gra2PXYdyx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dyx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ %...................................................................................
% $$$ figure(40)
% $$$ subplot(2,3,1)
% $$$ imagesc(xaxis,yaxis,dPXYdx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y) / dx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,2)
% $$$ imagesc(xaxis,yaxis,d2PXYdxx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dxx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,3)
% $$$ imagesc(xaxis,yaxis,d2PXYdxy)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dxy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,4)
% $$$ imagesc(xaxis,yaxis,dPXYdy)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y) / dy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,5)
% $$$ imagesc(xaxis,yaxis,d2PXYdyy)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dyy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,6)
% $$$ imagesc(xaxis,yaxis,d2PXYdyx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y) / dyx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ %...................................................................................
% $$$ %===================================================================================
% $$$ %...................................................................................
% $$$ figure(200)
% $$$ subplot(2,3,1)
% $$$ imagesc(xaxis,yaxis,delPXYalfadx)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y)^{\alpha} / dx')
% $$$ grid on
% $$$ subplot(2,3,2)
% $$$ imagesc(xaxis,yaxis,del2PXYalfadxx)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dxx')
% $$$ grid on
% $$$ subplot(2,3,3)
% $$$ imagesc(xaxis,yaxis,del2PXYalfadxy)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dxy')
% $$$ grid on
% $$$ subplot(2,3,4)
% $$$ imagesc(xaxis,yaxis,delPXYalfady)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y)^{\alpha} / dy')
% $$$ grid on
% $$$ subplot(2,3,5)
% $$$ imagesc(xaxis,yaxis,del2PXYalfadyy)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dyy')
% $$$ grid on
% $$$ subplot(2,3,6)
% $$$ imagesc(xaxis,yaxis,del2PXYalfadyx)
% $$$ colorbar('horizontal')
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dyx')
% $$$ grid on
% $$$ %...................................................................................
% $$$ figure(300)
% $$$ subplot(2,3,1)
% $$$ imagesc(xaxis,yaxis,graPXYalfadx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y)^{\alpha} / dx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,2)
% $$$ imagesc(xaxis,yaxis,gra2PXYalfadxx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dxx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,3)
% $$$ imagesc(xaxis,yaxis,gra2PXYalfadxy)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dxy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,4)
% $$$ imagesc(xaxis,yaxis,graPXYalfady)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y)^{\alpha} / dy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,5)
% $$$ imagesc(xaxis,yaxis,gra2PXYalfadyy)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dyy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,6)
% $$$ imagesc(xaxis,yaxis,gra2PXYalfadyx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dyx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ %...................................................................................
% $$$ figure(400)
% $$$ subplot(2,3,1)
% $$$ imagesc(xaxis,yaxis,dPXYalfadx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y)^{\alpha} / dx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,2)
% $$$ imagesc(xaxis,yaxis,d2PXYalfadxx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dxx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,3)
% $$$ imagesc(xaxis,yaxis,d2PXYalfadxy)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dxy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,4)
% $$$ imagesc(xaxis,yaxis,dPXYalfady)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d p(x,y)^{\alpha} / dy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,5)
% $$$ imagesc(xaxis,yaxis,d2PXYalfadyy)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dyy')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ subplot(2,3,6)
% $$$ imagesc(xaxis,yaxis,d2PXYalfadyx)
% $$$ xlabel('phy cell size (x)')
% $$$ ylabel('phy peak sst  (y)')
% $$$ title('d2 p(x,y)^{\alpha} / dyx')
% $$$ colorbar('horizontal')
% $$$ grid on
% $$$ %...................................................................................
% $$$ %===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GRAZING FUNCTIONAL RESPONSE KTW 2D:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%GENERAL GRAZING FUNCTIONAL RESPONSE KTW 2D:
%...................................................................................
Z = 1.0;
kgz = 0.5; 
gzmax = 1.0;
Vmax = gzmax*Z;
%...................................................................................
Qsat001 = (sumPXYdxdy / (sumPXYdxdy + kgz)); 
Qsat002 = (intPxydxdy / (intPxydxdy + kgz)); 
%...................................................................................
GXY001 = (PXYalfa ./ sumPXYalfadxdy) * Qsat001 * Vmax; %Numerical.
GXY002 = (PXYalfa ./ intPxyalfadxdy) * Qsat002 * Vmax; %Analytical.
%...................................................................................
GXY = GXY002; 
Qsat = Qsat002;
%...................................................................................
%===================================================================================
%BIOMASS SPECIFIC [d-1] GRAZING FUNCTIONAL RESPONSE KTW 2D:
%...................................................................................
gXY001 = GXY./PXY; 
gXY002 = FXY.^(alfa-1) * ((Qsat*Vmax)/(Ctot*Acte));
%...................................................................................
gXY = gXY002;
%...................................................................................
figure(100)
himg = pcolor(xaxis,yaxis,gXY);
set(himg,'edgecolor','none')
hbar = colorbar('vertical');
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('Grazing KTW -- g(x,y)')
%%return
%...................................................................................
%===================================================================================
%NUMERICAL DERIVATIVES OF BIOMASS SPECIFIC GRAZING FUNCTIONAL RESPONSE KTW 2D:
%...................................................................................
delgXYdx = diff(gXY,nth1,ndimX)/dx;
delgXYdy = diff(gXY,nth1,ndimY)/dy;
%...................................................................................
delgXYdx = [delgXYdx,ones(nsize,1)*nan]; 
delgXYdy = [delgXYdy;ones(1,msize)*nan]; 
%...................................................................................
del2gXYdxx = diff(delgXYdx,nth1,ndimX)/dx;
del2gXYdyy = diff(delgXYdy,nth1,ndimY)/dy;
%...................................................................................
del2gXYdxy = diff(delgXYdy,nth1,ndimX)/dx;
del2gXYdyx = diff(delgXYdx,nth1,ndimY)/dy;
%...................................................................................
[gragXYdx  ,gragXYdy]   = gradient(gXY     ,dx,dy);
[gra2gXYdxx,gra2gXYdyx] = gradient(gragXYdx,dx,dy);
[gra2gXYdxy,gra2gXYdyy] = gradient(gragXYdy,dx,dy);
%...................................................................................
%===================================================================================
%ANALYTICAL DERIVATIVES OF BIOMASS SPECIFIC GRAZING FUNCTIONAL RESPONSE KTW 2D:
%...................................................................................
Dcte = (Qsat*Vmax)/(Ctot*Acte);
%...................................................................................
dgXYdx000 = ( ((dFXYalfadx .* FXY) - (FXYalfa .* dFXYdx)) ./ (FXY.*FXY) ) * Dcte; %OKAY
dgXYdy000 = ( ((dFXYalfady .* FXY) - (FXYalfa .* dFXYdy)) ./ (FXY.*FXY) ) * Dcte; %OKAY
%...................................................................................
dgXYdx001 = ( FXY.^(alfa-1) .* SX * (alfa - 1) ) * Dcte; %OKAY
dgXYdy001 = ( FXY.^(alfa-1) .* SY * (alfa - 1) ) * Dcte; %OKAY
%...................................................................................
dgXYdx = dgXYdx000;
dgXYdy = dgXYdy000;
%...................................................................................
% $$$ dgXYdx = dgXYdx001;
% $$$ dgXYdy = dgXYdy001;
%...................................................................................
%===================================================================================
%...................................................................................
d2gXYdxx001 = ( (alfa-1) * FXY.^(alfa-1) .* ((alfa-1) * (SX).^2 - (1/sigmax^2)) ) * Dcte; %OKAY
d2gXYdyy001 = ( (alfa-1) * FXY.^(alfa-1) .* ((alfa-1) * (SY).^2 - (1/sigmay^2)) ) * Dcte; %OKAY
%...................................................................................
d2gXYdxx = d2gXYdxx001;
d2gXYdyy = d2gXYdyy001;
%...................................................................................
% $$$ figure(1)
% $$$ subplot(2,2,1)
% $$$ plot(dgXYdx(:),gragXYdx(:),'*')
% $$$ subplot(2,2,2)
% $$$ plot(dgXYdy(:),gragXYdy(:),'*')
% $$$ subplot(2,2,3)
% $$$ plot(d2gXYdxx(:),gra2gXYdxx(:),'*')
% $$$ subplot(2,2,4)
% $$$ plot(d2gXYdyy(:),gra2gXYdyy(:),'*')
%...................................................................................
%===================================================================================
%...................................................................................
d2gXYdxy = FXY.^(alfa-1) .* (SX.*SY) .* (alfa-1).*(alfa-1) .* Dcte;
d2gXYdyx = FXY.^(alfa-1) .* (SX.*SY) .* (alfa-1).*(alfa-1) .* Dcte;
%...................................................................................
%===================================================================================
%...................................................................................
figure(2000)
subplot(2,3,1)
imagesc(xaxis,yaxis,delgXYdx)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d g(x,y) / dx')
grid on
subplot(2,3,2)
imagesc(xaxis,yaxis,del2gXYdxx)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dxx')
grid on
subplot(2,3,3)
imagesc(xaxis,yaxis,del2gXYdxy)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dxy')
grid on
subplot(2,3,4)
imagesc(xaxis,yaxis,delgXYdy)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d g(x,y) / dy')
grid on
subplot(2,3,5)
imagesc(xaxis,yaxis,del2gXYdyy)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dyy')
grid on
subplot(2,3,6)
imagesc(xaxis,yaxis,del2gXYdyx)
colorbar('horizontal')
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dyx')
grid on
%...................................................................................
figure(3000)
subplot(2,3,1)
imagesc(xaxis,yaxis,gragXYdx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d g(x,y) / dx')
colorbar('horizontal')
grid on
subplot(2,3,2)
imagesc(xaxis,yaxis,gra2gXYdxx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dxx')
colorbar('horizontal')
grid on
subplot(2,3,3)
imagesc(xaxis,yaxis,gra2gXYdxy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dxy')
colorbar('horizontal')
grid on
subplot(2,3,4)
imagesc(xaxis,yaxis,gragXYdy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d g(x,y) / dy')
colorbar('horizontal')
grid on
subplot(2,3,5)
imagesc(xaxis,yaxis,gra2gXYdyy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dyy')
colorbar('horizontal')
grid on
subplot(2,3,6)
imagesc(xaxis,yaxis,gra2gXYdyx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dyx')
colorbar('horizontal')
grid on
%...................................................................................
figure(4000)
subplot(2,3,1)
imagesc(xaxis,yaxis,dgXYdx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d g(x,y) / dx')
colorbar('horizontal')
grid on
subplot(2,3,2)
imagesc(xaxis,yaxis,d2gXYdxx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dxx')
colorbar('horizontal')
grid on
subplot(2,3,3)
imagesc(xaxis,yaxis,d2gXYdxy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dxy')
colorbar('horizontal')
grid on
subplot(2,3,4)
imagesc(xaxis,yaxis,dgXYdy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d g(x,y) / dy')
colorbar('horizontal')
grid on
subplot(2,3,5)
imagesc(xaxis,yaxis,d2gXYdyy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dyy')
colorbar('horizontal')
grid on
subplot(2,3,6)
imagesc(xaxis,yaxis,d2gXYdyx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 g(x,y) / dyx')
colorbar('horizontal')
grid on
%...................................................................................
%===================================================================================
