%%function [] = jamstecrest_sstniches()

%===================================================================================
%...................................................................................
sst = [10:0.2:40]; %External sst forcing [degress C]
%...................................................................................
ymin = 0; %log([degrees C])
ymax = 40; %log([degrees C])
ydel = (ymax - ymin)/256; 
yaxis  = [ymin:ydel:ymax]; %Internal sst phy peak (maximum growth) [degress C]
%...................................................................................
yopt = sst; %Optimal value of internal phy sst peak should match the external sst forcing.
%...................................................................................
%%Gamma = (ymax - ymin)/6;
Gamma = (ymax - ymin)/4;
Theta = 2.0; %It *must* be [2, 4, 6, etc]
%...................................................................................
[YAXIS,YOPT] = meshgrid(yaxis,yopt);
[YAXIS,SST ] = meshgrid(yaxis,sst );
%...................................................................................
FY001 = (1.0) * exp(-(((YAXIS - YOPT).^Theta)/(2*Gamma^Theta))); %Pseudo-Gaussian curve reaching max.value = 1.0 (ie. area bigger than 1.0)
FY002 = (1.0) * exp(-(((YAXIS - SST ).^Theta)/(2*Gamma^Theta))); %Pseudo-Gaussian curve reaching max.value = 1.0 (ie. area bigger than 1.0)
%...................................................................................
FY = FY002; 
%...................................................................................
figure(10)
subplot(2,2,1)
surf(SST,YAXIS,FY) 
xlabel('External SST')
ylabel('Internal phy sst peak')
% $$$ set(gca,'Xdir','normal')
% $$$ set(gca,'Ydir','reverse')
shading interp 
grid on
subplot(2,2,2)
contourf(SST,YAXIS,FY)
xlabel('External SST')
ylabel('Internal phy sst peak')
grid on
%...................................................................................
%===================================================================================
%...................................................................................
xmin = log(1.0);
xmax = log(200.00);
xdel = ((xmax - xmin)/256); 
xaxis  = [xmin:xdel:xmax]; %[log(um)] 
%...................................................................................
din = [0:0.05:10.0];
%...................................................................................
[XAXIS,DIN] = meshgrid(xaxis,din);
%...................................................................................
mp = 0.20; %Phy background mortality [d-1]
mu0 = 1.0; %Phy maximum grazing rate [d-1] at log(ESD) = 1.0 
ks0 = 0.1; %Phy half-sat uptake [mmolN*m-3] at log(ESD) = 1.0 
amu = 1.00; %Size scaling factor for muPhy [n.d.]
aks = 2.00; %Size scaling factor for ksPhy [n.d.]
%...................................................................................
MUX001 = mu0 * exp(amu*XAXIS);
KSX001 = ks0 * exp(aks*XAXIS);
%...................................................................................
mumax = 3.0;
mux = mumax * (1 - exp(-0.02*exp(xaxis)));
MUX002 = mumax * (1 - exp(-0.02*exp(XAXIS)));
knp1 = 0.5; 
mup1 = 2.0; %[d-1]
alfa1 = mup1/knp1; %[m3*mmolN-1*d-1]
Cte = alfa1*mup1; %[m3*mmolN-1*d-2]
ksx = (mux.^2)/Cte;
KSX002 = (MUX002.^2)/Cte;
%...................................................................................
% $$$ MUX = MUX001; 
% $$$ KSX = KSX001; 
%...................................................................................
MUX = MUX002; 
KSX = KSX002; 
%...................................................................................
LX = (KSX ./ (KSX + DIN));
QX = (DIN ./ (DIN + KSX));
UX = MUX .* QX; 
%...................................................................................
figure(20)
subplot(2,2,1)
surf(DIN,XAXIS,UX) 
xlabel('DIN')
ylabel('phy cell size')
% $$$ set(gca,'Xdir','normal')
% $$$ set(gca,'Ydir','reverse')
shading interp 
grid on
subplot(2,2,2)
contourf(DIN,XAXIS,UX)
xlabel('DIN')
ylabel('phy cell size')
grid on
%...................................................................................
%===================================================================================
