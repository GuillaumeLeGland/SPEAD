function [ ] = jamstecrest_derivatives2D_uptakegrowth()
global mup0 knp0 
global amup aknp 

%===================================================================================
%CELL SIZE:
%...................................................................................
xmin = log(0.02); %USAR ESTOS!!!!
xmax = log(200.0);
xdel = (xmax - xmin)/(256-1); 
xaxis = [xmin:xdel:xmax]; 
%...................................................................................
din = [0.1:0.1:10.0]; 
%...................................................................................
%===================================================================================
%TEMPERATURE:
%...................................................................................
ymin = 0; %log([degrees C])
ymax = 50; %log([degrees C])
ydel = (ymax - ymin)/(128-1); 
yaxis  = [ymin:ydel:ymax]; %Internal sst phy peak (maximum growth) [degress C]
%...................................................................................
Gammay = (ymax - ymin)/6;
%...................................................................................
% $$$ Q10 = 1.0; %eppley-NOT
Q10 = 2.0; %eppley-YES
%...................................................................................
sst = [10:5:30]; %External sst forcing [degress C]
%...................................................................................
sst0 = 15.0;
%...................................................................................
%===================================================================================
%...................................................................................
[XAXIS,YAXIS] = meshgrid(xaxis,yaxis);
%...................................................................................
[mux,ksx] = jamstecrest_uptaketradeoff(xaxis,'Lanimal');
%...................................................................................
% $$$ lx = ones(length(xaxis),length(yaxis))*nan;
% $$$ qx = ones(length(xaxis),length(yaxis))*nan;
% $$$ ux = ones(length(xaxis),length(yaxis))*nan;
%...................................................................................
%===================================================================================
%UPTAKE RATE U(x,y) DISTRIBUTION 2D:
%...................................................................................
ssti = 25.0;
dinj =  0.5;
%...................................................................................
lx = (ksx ./ (ksx + dinj));
qx = (dinj ./ (dinj + ksx));
ux = qx .* mux; 
%...................................................................................
ly = - (yaxis - ssti) / Gammay^2; 
qy = exp( - (yaxis - ssti).^2 / (2*Gammay.^2) );
uy = qy * Q10^((ssti - sst0)/10); 
%...................................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% $$$ %...................................................................................
% $$$ xm = log(2.0);
% $$$ ym = ssti;
% $$$ %...................................................................................
% $$$ Gammax = (xmax - xmin) / 6; 
% $$$ Gammay = (ymax - ymin) / 6; 
% $$$ %...................................................................................
% $$$ lx = - (xaxis - xm) / Gammax^2; 
% $$$ ly = - (yaxis - ym) / Gammay^2; 
% $$$ %...................................................................................
% $$$ qx = exp( - (xaxis - xm).^2 / (2*Gammax.^2) );
% $$$ qy = exp( - (yaxis - ym).^2 / (2*Gammay.^2) );
% $$$ %...................................................................................
% $$$ % $$$ ux = qx .* (mup0); 
% $$$ % $$$ uy = qy * Q10^((ssti - sst0)/10); 
% $$$ %...................................................................................
% $$$ ux = qx;
% $$$ uy = qy;
% $$$ %...................................................................................
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%...................................................................................
[LX, ~] = meshgrid(lx,yaxis);
[QX, ~] = meshgrid(qx,yaxis);
%...................................................................................
[junkarg, LY] = meshgrid(xaxis,ly);
[junkarg, QY] = meshgrid(xaxis,qy);
%...................................................................................
[UX,UY] = meshgrid(ux,uy);
%...................................................................................
UXY = UX .* UY; 
%...................................................................................
figure(10)
% $$$ contourf(XAXIS,YAXIS,UXY);
% $$$ himg = imagesc(UXY,[0 10.0]);
himg = pcolor(XAXIS,YAXIS,UXY);
%%himg = pcolor(XAXIS,YAXIS,log2(UXY));
set(himg,'edgecolor','none')
%%caxis([log2(0.01) log2(15.0)]) %put before colorbar
hbar = colorbar('vertical');
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('Uptake -- u(x,y)')
%%set(hbar,'Ylim',[0 10])
%...................................................................................
%===================================================================================
%NUMERICAL DERIVATIVES OF U(x,y):
%...................................................................................
nth1 = 1; %First derivative (gradient)
nth2 = 2; %Second derivative (laplacian)
%...................................................................................
ndimY = 1; %diff between rows (ie. gradient in y-direction or columnwise)
ndimX = 2; %diff between cols (ie. gradient in x-direction or rowes-wise)
%...................................................................................
delUdx   = diff(UXY,nth1,ndimX)   /xdel;
delUdy   = diff(UXY,nth1,ndimY)   /ydel;
del2Udxx = diff(UXY,nth2,ndimX)   /xdel;
del2Udyy = diff(UXY,nth2,ndimY)   /ydel;
del2Udxy = diff(delUdy,nth1,ndimX)/xdel;
del2Udyx = diff(delUdx,nth1,ndimY)/ydel;
%...................................................................................
[graUdx,graUdy]     = gradient(UXY   ,xdel,ydel);
[gra2Udxx,gra2Udyx] = gradient(graUdx,xdel,ydel);
[gra2Udxy,gra2Udyy] = gradient(graUdy,xdel,ydel);
%...................................................................................
%===================================================================================
%ANALYTICAL DERIVATIVES OF U(x,y):
%...................................................................................
dUdx = UXY .* (amup - aknp * LX); 
dUdy = UXY .* LY; 
%...................................................................................
d2Udxx = UXY .* ( (amup - aknp * LX).^2 - (aknp^2 * LX .* QX) );
d2Udyy = UXY .* ( (LY).^2 - (1/Gammay^2) );
%...................................................................................
d2Udxy = UXY .* (amup - aknp * LX) .* LY; 
d2Udyx = d2Udxy; 
%...................................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%...................................................................................
% $$$ dUdx = UXY .* LX; 
% $$$ dUdy = UXY .* LY; 
%...................................................................................
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%...................................................................................
% $$$ figure(50)
% $$$ subplot(2,2,1)
% $$$ plot(delUdx(:,:),dUdx(:,1:end-1),'*')
% $$$ set(gca,'Xlim',[-1.00 1.00],'Ylim',[-1.00 1.00])
% $$$ grid on
% $$$ subplot(2,2,2)
% $$$ plot(delUdy(:,:),dUdy(1:end-1,:),'*')
% $$$ set(gca,'Xlim',[-0.10 0.10],'Ylim',[-0.10 0.10])
% $$$ grid on
% $$$ subplot(2,2,3)
% $$$ plot(graUdx(:,:),dUdx,'*')
% $$$ set(gca,'Xlim',[-1.00 1.00],'Ylim',[-1.00 1.00])
% $$$ grid on
% $$$ subplot(2,2,4)
% $$$ plot(graUdy(:,:),dUdy,'*')
% $$$ set(gca,'Xlim',[-0.10 0.10],'Ylim',[-0.10 0.10])
% $$$ grid on
% $$$ return
%...................................................................................
%===================================================================================
%PLOTS:
figure(20)
subplot(2,3,1)
imagesc(xaxis,yaxis,delUdx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d u(x,y) / dx')
colorbar('horizontal')
grid on
subplot(2,3,2)
imagesc(xaxis,yaxis,del2Udxx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dxx')
colorbar('horizontal')
grid on
subplot(2,3,3)
imagesc(xaxis,yaxis,del2Udxy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dxy')
colorbar('horizontal')
grid on
subplot(2,3,4)
imagesc(xaxis,yaxis,delUdy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d u(x,y) / dy')
colorbar('horizontal')
grid on
subplot(2,3,5)
imagesc(xaxis,yaxis,del2Udyy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dyy')
colorbar('horizontal')
grid on
subplot(2,3,6)
imagesc(xaxis,yaxis,del2Udyx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dyx')
colorbar('horizontal')
grid on
%...................................................................................
figure(30)
subplot(2,3,1)
imagesc(xaxis,yaxis,graUdx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d u(x,y) / dx')
colorbar('horizontal')
grid on
subplot(2,3,2)
imagesc(xaxis,yaxis,gra2Udxx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dxx')
colorbar('horizontal')
grid on
subplot(2,3,3)
imagesc(xaxis,yaxis,gra2Udxy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dxy')
colorbar('horizontal')
grid on
subplot(2,3,4)
imagesc(xaxis,yaxis,graUdy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d u(x,y) / dy')
colorbar('horizontal')
grid on
subplot(2,3,5)
imagesc(xaxis,yaxis,gra2Udyy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dyy')
colorbar('horizontal')
grid on
subplot(2,3,6)
imagesc(xaxis,yaxis,gra2Udyx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dyx')
colorbar('horizontal')
grid on
%...................................................................................
figure(40)
subplot(2,3,1)
imagesc(xaxis,yaxis,dUdx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d u(x,y) / dx')
colorbar('horizontal')
grid on
subplot(2,3,2)
imagesc(xaxis,yaxis,d2Udxx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dxx')
colorbar('horizontal')
grid on
subplot(2,3,3)
imagesc(xaxis,yaxis,d2Udxy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dxy')
colorbar('horizontal')
grid on
subplot(2,3,4)
imagesc(xaxis,yaxis,dUdy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d u(x,y) / dy')
colorbar('horizontal')
grid on
subplot(2,3,5)
imagesc(xaxis,yaxis,d2Udyy)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dyy')
colorbar('horizontal')
grid on
subplot(2,3,6)
imagesc(xaxis,yaxis,d2Udyx)
xlabel('phy cell size (x)')
ylabel('phy peak sst  (y)')
title('d2 u(x,y) / dyx')
colorbar('horizontal')
grid on
%...................................................................................
% $$$ return
%===================================================================================
% $$$ for i = 1:length(sst)
% $$$     for j = 1:length(din)
% $$$ 	%...........................................................................
% $$$ 	ssti = sst(i);
% $$$ 	dinj = din(j);
% $$$ 	%...........................................................................
% $$$ 	lx = (ksx ./ (ksx + dinj));
% $$$ 	qx = (dinj ./ (dinj + ksx));
% $$$ 	ux = qx .* mux; 
% $$$ 	%...........................................................................
% $$$ 	qy = exp( - (yaxis - ssti).^2 / (2*Gammay.^2) );
% $$$ 	uy = qy * Q10^((ssti - sst0)/10); 
% $$$ 	%...........................................................................
% $$$ 	[UX,UY] = meshgrid(ux,uy);
% $$$ 	%...........................................................................
% $$$ 	UXY = UX .* UY; 
% $$$ 	%...........................................................................
% $$$ 	figure(20)
% $$$ 	%...........................................................................
% $$$ % $$$ 	contourf(XAXIS,YAXIS,UXY)
% $$$ % $$$ 	himg = imagesc(UXY,[0 10.0])
% $$$ 	himg = pcolor(XAXIS,YAXIS,log2(UXY));
% $$$ 	set(himg,'edgecolor','none')
% $$$ 	caxis([log2(0.01) log2(15.0)]) %put before colorbar
% $$$ 	hbar = colorbar('vertic');
% $$$ 	xlabel('phy cell size (x)')
% $$$ 	ylabel('phy peak sst  (y)')
% $$$ 	title('Uptake -- u(x,y)')
% $$$ % $$$ 	set(hbar,'Ylim',[0 10])
% $$$ % $$$ 	pause
% $$$ 	%...........................................................................
% $$$     end
% $$$ end
%===================================================================================
