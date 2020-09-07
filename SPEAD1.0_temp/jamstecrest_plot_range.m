% Script to plot the minimum, average and maximum of all state variables, for 
% different trait diffusivity values (Le Gland, 19/07/2020)

%function [hfig] = jamstecrest_barplot_range (P,xave,yave,xstd,ystd,xycor,fignum,mypackages)

% Different trait diffusivity values
nux = [0, 0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1];

P = [146.46, 146.94, 146.85, 146.81, 146.84, 146.9, 146.62, 144.89, 140.62, 130.8];

xave_min = [ -0.56398, -0.48275, -0.47482, -0.51174, -0.60947, -0.82977, -1.0839, -1.2877, -1.3685, -1.2668];
xave_ave = [ -0.1434, -0.044702, -0.046165, -0.049356, -0.05576, -0.065778, -0.068322, -0.050009, -0.008878, 0.067136];
xave_max = [ 0.18996, 0.32952, 0.3502, 0.39681, 0.46958, 0.59707, 0.74928, 0.90348, 0.99617, 1.0434];

yave_min = [ 21.77, 22.125, 22.166, 22.193, 22.153, 21.982, 21.722, 21.425, 21.219, 21.06];
yave_ave = [ 23.097, 23.397, 23.376, 23.356, 23.321, 23.264, 23.204, 23.144, 23.107, 23.099];
yave_max = [ 24.774, 24.977, 25.098, 25.25, 25.526, 26.085, 26.685, 27.864, 28.593, 29.048];

xstd_min = [0.11902, 0.14004, 0.16653, 0.19949, 0.24262, 0.30846, 0.39435, 0.52404, 0.67869, 0.89226];
xstd_ave = [0.19703, 0.2179, 0.23232, 0.2591, 0.29797, 0.36117, 0.43998, 0.55639, 0.69936, 0.90531];
xstd_max = [0.41493, 0.45569, 0.46993, 0.52211, 0.61532, 0.77059, 0.91397, 1.0296, 1.1096, 1.1884];

ystd_min = [0.47421, 0.47997, 0.53662, 0.59346, 0.67275, 0.81013, 1.0199, 1.4053, 1.934, 2.6882];
ystd_ave = [0.78511, 0.76384, 0.82253, 0.88001, 0.9737, 1.1526, 1.4088, 1.8078, 2.2984, 2.9974];
ystd_max = [1.6528, 1.5819, 1.5899, 1.6186, 1.7176, 1.921, 2.149, 2.4197, 2.7342, 3.2874];

xycor_min = [-0.999, -0.98496, -0.92442, -0.86919, -0.82218, -0.77659, -0.7115, -0.58179, -0.41653, -0.24811];
xycor_ave = [-0.99859, -0.92645, -0.67944, -0.49249, -0.34291, -0.21532, -0.13009, -0.064981, -0.02881, -0.0095926];
xycor_max = [-0.99426, -0.85802, -0.44918, -0.21148, -0.066757, -0.0090818, -0.00029523, 0.028865, 0.0050423, 0.0032445];

hfig = figure(90);

hplot = subplot(2,3,1);
bar(P,'g')
%plot(P,'g.','markersize',25)
xlabel(hplot,'Trait diffusivity x [n. d.]')
ylabel(hplot,'Primary production [gC.m^{-2}.yr^{-1}]')
set(hplot,'Ylim',[130 150],'Xlim',[0 11],'Xtick',[1:1:10],'XtickLabel',nux)
grid on

hplot = subplot(2,3,2);
Ytl = [0.082, 0.14, 0.22, 0.37, 0.61, 1, 1.65, 2.72, 4.48]; % Half-saturation labels 
plot(xave_ave,'b.','markersize',30)
hold on
plot([0 11],[1.2 1.2],'k-','linewidth',3)
plot([0 11],[-1.7 -1.7],'k-','linewidth',3)
xlabel(hplot,'Trait diffusivity x [n. d.]')
ylabel(hplot,'Mean half-saturation [mmolN.m^{-3}]')
set(hplot,'Ylim',[-2.5 1.5],'Xlim',[0 11],'Ytick',[-2.5:0.5:1.5],'YtickLabel',Ytl,'Xtick',[1:1:10],'XtickLabel',nux)
ha = errorbar([1:1:10],xave_ave,xave_ave-xave_min,xave_max-xave_ave,'b.','linewidth',2);
% To change the width of caps (in R2016b and forward, errorbar has a property called 'capsize')
hb = get(ha,'children'); 
Xdata = get(hb(2),'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
% xleft and xright contain the indices of the left and right
%  endpoints of the horizontal lines
xleft = temp; xright = temp+1;
% Increase line length by 0.2 units
Xdata(xleft) = Xdata(xleft) - .15;
Xdata(xright) = Xdata(xright) + .15;
set(hb(2),'Xdata',Xdata)
grid on

hplot = subplot(2,3,3); 
plot(yave_ave,'r.','markersize',30)
hold on
plot([0 11],[27.8 27.8],'k-','linewidth',3)
plot([0 11],[18.5 18.5],'k-','linewidth',3)
xlabel(hplot,'Trait diffusivity x [n. d.]')
ylabel(hplot,'Optimal temperature [°C]')
set(hplot,'Ylim',[18 30],'Xlim',[0 11],'Ytick',[18:2:30],'Xtick',[1:1:10],'XtickLabel',nux)
% To change the width of caps (in R2016b and forward, errorbar has a property called 'capsize')
ha = errorbar([1:1:10],yave_ave,yave_ave-yave_min,yave_max-yave_ave,'r.','linewidth',2);
hb = get(ha,'children'); 
Xdata = get(hb(2),'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
% xleft and xright contain the indices of the left and right
%  endpoints of the horizontal lines
xleft = temp; xright = temp+1;
% Increase line length by 0.2 units
Xdata(xleft) = Xdata(xleft) - .15;
Xdata(xright) = Xdata(xright) + .15;
set(hb(2),'Xdata',Xdata)
grid on

hplot = subplot(2,3,5); 
plot(xstd_ave,'b.','markersize',30)
hold on
plot([0 11],[2/sqrt(3) 2/sqrt(3)],'k-','linewidth',3)
xlabel(hplot,'Trait diffusivity x [n. d.]')
ylabel(hplot,'Standard deviation of logKn')
set(hplot,'Ylim',[0 1.4],'Xlim',[0 11],'Ytick',[0:0.2:1.4],'Xtick',[1:1:10],'XtickLabel',nux)
ha = errorbar([1:1:10],xstd_ave,xstd_ave-xstd_min,xstd_max-xstd_ave,'b.','linewidth',2);
% To change the width of caps (in R2016b and forward, errorbar has a property called 'capsize')
hb = get(ha,'children'); 
Xdata = get(hb(2),'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
% xleft and xright contain the indices of the left and right
%  endpoints of the horizontal lines
xleft = temp; xright = temp+1;
% Increase line length by 0.2 units
Xdata(xleft) = Xdata(xleft) - .15;
Xdata(xright) = Xdata(xright) + .15;
set(hb(2),'Xdata',Xdata)
grid on

hplot = subplot(2,3,6); 
plot(ystd_ave,'r.','markersize',30)
hold on
plot([0 11],[6/sqrt(3) 6/sqrt(3)],'k-','linewidth',3)
xlabel(hplot,'Trait diffusivity x [n. d.]')
ylabel(hplot,'Standard deviation of optimal temperature')
set(hplot,'Ylim',[0 4],'Xlim',[0 11],'Ytick',[0:0.5:4],'Xtick',[1:1:10],'XtickLabel',nux)
ha = errorbar([1:1:10],ystd_ave,ystd_ave-ystd_min,ystd_max-ystd_ave,'r.','linewidth',2);
% To change the width of caps (in R2016b and forward, errorbar has a property called 'capsize')
hb = get(ha,'children'); 
Xdata = get(hb(2),'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
% xleft and xright contain the indices of the left and right
%  endpoints of the horizontal lines
xleft = temp; xright = temp+1;
% Increase line length by 0.2 units
Xdata(xleft) = Xdata(xleft) - .15;
Xdata(xright) = Xdata(xright) + .15;
set(hb(2),'Xdata',Xdata)
grid on

hplot = subplot(2,3,4); 
plot(xycor_ave,'k.','markersize',30)
hold on
plot([0 11],[0 0],'k-','linewidth',3)
xlabel(hplot,'Trait diffusivity x [n. d.]')
ylabel(hplot,'Inter-trait covariance')
set(hplot,'Ylim',[-1 1],'Xlim',[0 11],'Ytick',[-1:0.25:1],'Xtick',[1:1:10],'XtickLabel',nux)
ha = errorbar([1:1:10],xycor_ave,xycor_ave-xycor_min,xycor_max-xycor_ave,'k.','linewidth',2);
% To change the width of caps (in R2016b and forward, errorbar has a property called 'capsize')
hb = get(ha,'children'); 
Xdata = get(hb(2),'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
% xleft and xright contain the indices of the left and right
%  endpoints of the horizontal lines
xleft = temp; xright = temp+1;
% Increase line length by 0.2 units
Xdata(xleft) = Xdata(xleft) - .15;
Xdata(xright) = Xdata(xright) + .15;
set(hb(2),'Xdata',Xdata)
grid on

%end
