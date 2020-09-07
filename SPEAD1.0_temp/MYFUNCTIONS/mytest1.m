function []=mytest1(yinput,xinput)

%........................
% $$$ [Xstar,xz,xzbitmap,xzbitmapinf,xzbitmapsup]=mybitmapscale(xinput);
% $$$ [Ystar,yz,yzbitmap,yzbitmapinf,yzbitmapsup]=mybitmapscale(yinput);
%........................
% $$$ [Xstar,xz,xzbitmap,xzbitmapinf,xzbitmapsup]=mybitmapscale(xinput,0.1);
% $$$ [Ystar,yz,yzbitmap,yzbitmapinf,yzbitmapsup]=mybitmapscale(yinput,0.1);
%........................
[Xstar,xz,xzbitmap,xzbitmapinf,xzbitmapsup]=mybitmapscale(xinput,0.5);
[Ystar,yz,yzbitmap,yzbitmapinf,yzbitmapsup]=mybitmapscale(yinput,0.5);
%........................

log10x=Xstar;
log10y=Ystar;

display('=== pto1 ===')
xZandZbitmap=[xz;xzbitmap]
yZandZbitmap=[yz;yzbitmap]
display('=== pto2 ===')

%.........................
x=Xstar(:);
y=Ystar(:);
[ychap,a,b,bs,aint,bint,res,resint,r2,F,pval,nptos]=sregress(y,x);
I=find(isnan(x)==0);
xx=x(I);
xmin=min(x);
xmax=max(x);
dx=(xmax-xmin)/100;
xxi=[xmin:dx:xmax];
J=find(isnan(x)==0 & isnan(ychap)==0);
yychap=ychap(J);
xx=x(J);
ychapI=interp1(xx,yychap,xxi);
%.........................

figure(9000)
hp1=plot(log10x,log10y,'*');
hold on
hp2=plot(xxi,ychapI,'k-');
hold off
set(gca,'Xlim',[xzbitmap(1) xzbitmap(end)],'Xtick',xzbitmap(:),'Xticklabel',xz)
set(gca,'Ylim',[yzbitmap(1) yzbitmap(end)],'Ytick',yzbitmap(:),'Yticklabel',yz)
grid on
