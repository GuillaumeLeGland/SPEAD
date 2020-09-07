function []=myloglogplot(Yvarname,Xvarname,yinput,xinput,fignum) %ver "myRegressLogPlot.m"!!!!
global iimin jjmin iimax jjmax

%************************************************************************************
% myloglogplot(Yvarname,Xvarname,yinput,xinput,fignum) %ver "myRegressLogPlot.m"!!!!
%************************************************************************************
%..........
% $$$ [yinput(:),xinput(:)]
% $$$ pause(1)
%..........
nptos=length(xinput);

%............................
Ix=find(abs(xinput)>0); %Get rid off zeros.
Iy=find(abs(yinput)>0);
xx=ones(nptos,1)*nan;
yy=ones(nptos,1)*nan;
xx(Ix)=xinput(Ix);
yy(Iy)=yinput(Iy);
x=xx;
y=yy;
xn=x./abs(nanmean(x(:)));
yn=y./nanmean(y(:));
%............................
logxinput=log10(xinput);
logyinput=log10(yinput);
[r2,corte,pente,npuntos]=myregress(logyinput,logxinput)
display('................')
xmin=min(x);
xmax=max(x);
dx=(xmax-xmin)/100;
xinterp=[xmin:dx:xmax];

a=10^corte;
b=pente;
ychap=a*xinterp.^b;
%............................
[logXstar,xz,xzbitmap,xzbitmapinf,xzbitmapsup]=mybitmapscale(xinput); %logXstar[0-256]
[logYstar,yz,yzbitmap,yzbitmapinf,yzbitmapsup]=mybitmapscale(yinput);

nlogXstar=logXstar/mean(logXstar(:));
nlogYstar=logYstar/mean(logYstar(:));

% $$$ [xinput(:),yinput(:)]
% $$$ [logXstar(:),logYstar(:)]
% $$$ [nlogXstar(:),nlogYstar(:)]
% $$$ pause

%............................
    
