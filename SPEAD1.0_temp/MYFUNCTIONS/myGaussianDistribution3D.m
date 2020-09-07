function [FXYZ,XAXIS,YAXIS,ZAXIS] = myGaussianDistribution3D(xmin,ymin,zmin,xmax,ymax,zmax,xlag,ylag,zlag,xsigma,ysigma,zsigma,xnbins,ynbins,znbins)

%**********************************************************************
%Use: [GaussXY]=myGaussDistribution(xrng,yrng,sigmaPcnt,xlag,ylag);
%
%----------------------------------------------------------------------
%sigmaXci = sqrt(0.01*xmax);
%GaussX = (1/(sigmaX*sqrt(2*pi))) * exp(-(xrng-xmean).^2/(2*sigmaX^2));
%----------------------------------------------------------------------
%**********************************************************************

%===================================================================================
%...................................................................................
% $$$ keymeshcall = 'NDGRID';
% $$$ keymeshcall = 'MESHGRID';
%...................................................................................
%===================================================================================
%...................................................................................
dx = ((xmax-xmin)/(xnbins-1)); 
dy = ((ymax-ymin)/(ynbins-1)); 
dz = ((zmax-zmin)/(znbins-1)); 
%...................................................................................
xaxis = [xmin:dx:xmax]; %[log(um)] 
yaxis = [ymin:dy:ymax]; %[log(um)] 
zaxis = [zmin:dz:zmax]; %[log(um)] 
%...................................................................................
xmean = mean(xaxis) + xlag;
ymean = mean(yaxis) + ylag;
zmean = mean(zaxis) + zlag;
%...................................................................................
% $$$ [XAXIS,YAXIS] = meshgrid(xaxis,yaxis); %Needs to do the transposed orientation.
%...................................................................................
[XAXIS,YAXIS,ZAXIS] = ndgrid(xaxis,yaxis,zaxis); %Transposed orientation of MESHGRID.
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%
%GAUSSIAN SURFACE:
%%%%%%%%%%%%%%%%%%
%===================================================================================
%-----------------------------------------------------------------------------------
% fx = 1 / (2*pi*xsigma*ysigma) * exp(-( (XAXIS - xmean).^2 / (2*xsigma^2) + (YAXIS - ymean).^2 / (2*ysigma^2) ));
%-----------------------------------------------------------------------------------
% fx = 1 / (sqrt((2*pi)^ndim*det(covxy))) * exp(-( (XAXIS - xmean).^2 / (2*xsigma^2) + (YAXIS - ymean).^2 / (2*ysigma^2) ));
%-----------------------------------------------------------------------------------
ndim = 3; 
%...................................................................................
covxyz = zeros(ndim,ndim); 
%...................................................................................
covxyz(1,1) = xsigma.^2; %Autocorrelations.
covxyz(2,2) = ysigma.^2; 
covxyz(3,3) = zsigma.^2; 
%...................................................................................
% $$$ covxyz(1,2) = 0.1^2; %Crosscorrelations (keep uncommented for "zero" crosscorrelations).
% $$$ covxyz(2,1) = 0.1^2; 
%...................................................................................
% $$$ Ax = 1 / (xsigma*(sqrt(2*pi))); %OKAY.
% $$$ Ay = 1 / (ysigma*(sqrt(2*pi))); %OKAY.
%...................................................................................
Ax  = 1 / (sqrt(2*pi*xsigma*xsigma)); %OKAY.
Ay  = 1 / (sqrt(2*pi*ysigma*ysigma)); %OKAY.
Az  = 1 / (sqrt(2*pi*zsigma*zsigma)); %OKAY.
%...................................................................................
Axyz = 1 / (sqrt((2*pi)^ndim*det(covxyz))); %OKAY.
%...................................................................................
Bx = (XAXIS-xmean).^2 / (2*xsigma^2);
By = (YAXIS-ymean).^2 / (2*ysigma^2);
Bz = (ZAXIS-zmean).^2 / (2*zsigma^2);
%...................................................................................
fx = Ax * exp(-(xaxis-xmean).^2 / (2*xsigma^2));
fy = Ay * exp(-(yaxis-ymean).^2 / (2*ysigma^2));
fz = Az * exp(-(zaxis-zmean).^2 / (2*zsigma^2));
%...................................................................................
FX = Ax * exp(-Bx);
FY = Ay * exp(-By);
FZ = Az * exp(-Bz);
%...................................................................................
FXYZ001 = FX .* FY .* FZ; 
%...................................................................................
FXYZ002 = Axyz * exp(-(Bx + By + Bz)); 
%...................................................................................
FXYZ = FXYZ001;
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE THE INTEGRAL UNDER THE PDF SURFACE TO BE SURE THAT IS EQUAL TO ONE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
[sumFXYZdxdydz] = myintegralpdf3D(FXYZ,xmax,xmin,ymax,ymin,zmax,zmin,xnbins,ynbins,znbins);
%...................................................................................
dxdydz = dx*dy*dz;
%...................................................................................
FXYZdxdydz = FXYZ * dxdydz;
%...................................................................................
sumFXYZdxdydz002 = sum(FXYZdxdydz); 
%...................................................................................
%===================================================================================

%***********************************************************************************
return
%...................................................................................
figure(10)
subplot(2,2,1)
plot(xaxis,fx,'*')
title('fx - Nopt')
subplot(2,2,2)
plot(yaxis,fy,'*')
title('fy - Topt')
subplot(2,2,3)
plot(zaxis,fz,'*')
title('fz - Iopt')
%...................................................................................
figure(15)
subplot(2,2,1)
plot(FXYZ001(:),FXYZ002(:),'*')
%...................................................................................
figure(20)
S001 = FXYZ001(:)/max(FXYZ001(:))*100;
C001 = FXYZ001(:);
subplot(2,2,1)
hplot = scatter3(XAXIS(:),YAXIS(:),ZAXIS(:),S001,C001,'filled');
%...................................................................................
S002 = FXYZ002(:)/max(FXYZ002(:))*100;
C002 = FXYZ002(:);
subplot(2,2,1)
hplot = scatter3(XAXIS(:),YAXIS(:),ZAXIS(:),S002,C002,'filled');
%...................................................................................
myTitle{1,1} = 'PAR opt 12C';
myTitle{2,1} = 'PAR opt 16C';
myTitle{3,1} = 'PAR opt 18C';
myTitle{4,1} = 'PAR opt 24C';
%...................................................................................
figure(25)
for j = 1:4
    subplot(2,2,j)
    mypcolorimagesc(xaxis,yaxis,FXYZ001(:,:,j))
    set(gca,'Clim',[0 max(FXYZ001(:))])
    xlabel('DIN opt')
    ylabel('SST opt')
    title(myTitle(j))
    shading flat 
    mycolorbar
end
%...................................................................................
figure(30)
for j = 1:4
    subplot(2,2,j)
    mypcolorimagesc(xaxis,yaxis,FXYZ002(:,:,j))
    set(gca,'Clim',[0 max(FXYZ002(:))])
    xlabel('DIN opt')
    ylabel('SST opt')
    title(myTitle(j))
    shading flat 
    mycolorbar
end
%...................................................................................
figure(35)
for j = 1:4
    subplot(2,2,j)
    hbar = bar3(FXYZ(:,:,j));
    for k = 1:length(hbar)
	zdata = get(hbar(k),'ZData');
	set(hbar(k),'CData',zdata)
	set(hbar(k),'FaceColor','interp')
	%%set(hbar(k),'LineWidth',[0.5])
	%%set(hbar(k),'Linestyle','none')
	%%set(hbar(k),'EdgeColor',[1 1 1])
    end
    xlabel('DIN opt')
    ylabel('SST opt')
    title(myTitle(j))
end
%...................................................................................
%===================================================================================
%***********************************************************************************
return