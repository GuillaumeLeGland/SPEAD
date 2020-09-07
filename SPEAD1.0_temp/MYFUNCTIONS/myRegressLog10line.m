function [hp1,hp2]=myRegressLog10line(xinput,yinput,varargin)

%***************************************************    
%Function: MYREGRESSLOG10LINE.m
%
%Use: [hp1,hp2]=myRegressLog10line(xinput,yinput)
%Use: [hp1,hp2]=myRegressLog10line(xinput,yinput,[],'yes')
%Use: [hp1,hp2]=myRegressLog10line(xinput,yinput,Wxy,'yes')
%Use: [hp1,hp2]=myRegressLog10line(xinput,yinput,Wxy,'yes','Robust')
%
%Inputs:
%
%xinput      = Vector (or array) of X-variable.
%yinput      = Vector (or array) of Y-variable.
%Wxy(1,1)    = xmin.
%Wxy(1,2)    = xmax.
%Wxy(2,1)    = ymin.
%Wxy(2,1)    = ymax.
%RegressLine = 'yes/no'.
%RegressFit  = 'Robust/Normal'.
%
%Outpus:
%
%hp1 = Plot handle for "dots".
%hp2 = Plot handle for "line".
%***************************************************    
%===================================================
nvarargin=length(varargin);
%===================================================
% $$$ disp('======= new data =========')
% $$$ xinput
% $$$ yinput
%===================================================
I=find(xinput<0);
J=find(yinput<0);
keyXdataNegative='no';
keyYdataNegative='no';
if length(I)>0
    %...........................................
    keyXdataNegative='yes';
    %...........................................
    disp('Ojo! some xinput data are negative!!')
    disp('I will force them to be positive (check if thats ok)')
    disp('Otherwise I can not take log10 of those date')
    %...........................................
    xinput=abs(xinput);
    if nvarargin>=1
	%...................
	xyminmax=varargin{:,1}; %Wxy
	%...................
	xmax=xyminmax(1,1)*(-1); %Change min/max "x" order.
	xmin=xyminmax(1,2)*(-1);
	%...................
	ymin=xyminmax(2,1);
	ymax=xyminmax(2,2);
	%...................
	Wxy=[xmin,xmax;ymin,ymax];
	%...................
	varargin{:,1}=Wxy;
	%...................
    end
    %...........................................
end
if length(J)>0
    %...........................................
    keyYdataNegative='yes';
    %...........................................
    disp('Ojo! some yinput data are negative!!')
    disp('I will force them to be positive (check if thats ok)')
    disp('Otherwise I can not take log10 of those date')
    %...........................................
    yinput=abs(yinput);
    if nvarargin>=1
	%...................
	xyminmax=varargin{:,1}; %Wxy
	%...................
	xmin=xyminmax(1,1);
	xmax=xyminmax(1,2);
	%...................
	ymax=xyminmax(2,1)*(-1); %Change min/max "y" order.
	ymin=xyminmax(2,2)*(-1);
	%...................
	Wxy=[xmin,xmax;ymin,ymax];
	%...................
	varargin{:,1}=Wxy;
	%...................
    end
    %...........................................
end
%===================================================
% $$$ xinput
% $$$ yinput
% $$$ pause
%===================================================
if nvarargin==0
    %...................
    RegressLine='yes';
% $$$     RegressLine='no';
    %...................
    RegressFit='Normal';
% $$$     RegressFit='Robust';
    %...................
    xmin=min(xinput(:));
    xmax=max(xinput(:));
    ymin=min(yinput(:));
    ymax=max(yinput(:));
    %...................
elseif nvarargin==1
    %...................
    xyminmax=varargin{:,1} %Wxy
    %...................
    xmin=xyminmax(1,1);
    xmax=xyminmax(1,2);
    ymin=xyminmax(2,1);
    ymax=xyminmax(2,2);
    %...................
elseif nvarargin==2
    %...................
    xyminmax=varargin{:,1} %Wxy
    %...................
    if isempty(xyminmax)==1
	xmin=min(xinput(:));
	xmax=max(xinput(:));
	ymin=min(yinput(:));
	ymax=max(yinput(:));
    elseif isempty(xyminmax)==0
	xmin=xyminmax(1,1);
	xmax=xyminmax(1,2);
	ymin=xyminmax(2,1);
	ymax=xyminmax(2,2);
    end
    %...................
    RegressLine=varargin{:,2}
    RegressFit='Normal';
    %...................
elseif nvarargin==3
    %...................
    xyminmax=varargin{:,1} %Wxy
    %...................
    if isempty(xyminmax)==1
	xmin=min(xinput(:));
	xmax=max(xinput(:));
	ymin=min(yinput(:));
	ymax=max(yinput(:));
    elseif isempty(xyminmax)==0
	xmin=xyminmax(1,1);
	xmax=xyminmax(1,2);
	ymin=xyminmax(2,1);
	ymax=xyminmax(2,2);
    end
    %...................
    RegressLine=varargin{:,2}
    RegressFit=varargin{:,3}
    %...................
end
xminmax=[xmin,xmax]
yminmax=[ymin,ymax]
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO AVOID EQUAL VALUES AND ADD A "MINI" DELTA INCREASE/DECREASE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(xinput);
xrand=rand(m,n);
yrand=rand(m,n);
fc=sqrt(eps);
xxinput=xinput+(xrand*fc.*xinput);
yyinput=yinput+(yrand*fc.*yinput);

%?????????????????????????????????
% $$$ varx=[xinput(:),xxinput(:)],pause
% $$$ vary=[yinput(:),yyinput(:)],pause
% $$$ sortx=sort([xinput(:),xxinput(:)]),pause
% $$$ sorty=sort([yinput(:),yyinput(:)]),pause
%?????????????????????????????????

xinput=xxinput;
yinput=yyinput;

%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%OBTAIN LOG DATA AND TICKS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%..........................
[Xstar,xz,xzbitmap,xzbitmapinf,xzbitmapsup]=mybitmapscale(xinput,'Log10',xmin,xmax);
[Ystar,yz,yzbitmap,yzbitmapinf,yzbitmapsup]=mybitmapscale(yinput,'Log10',ymin,ymax);
%..........................
log10x=Xstar;
log10y=Ystar;
%..........................

%REGRESSION OVER LOG DATA:
if strcmp(RegressLine,'yes')
    %.........................
    x=Xstar(:);
    y=Ystar(:);
    %.........................
    if strcmp(RegressFit,'Normal')
	[a,b,bs,r2,ychap,nptos]=myRegress(y,x);
	xmin=min(x);
	xmax=max(x);
	dx=(xmax-xmin)/100;
	xchapI=[xmin:dx:xmax];
	ychapI = a + (b*xchapI);
	%.........................
    elseif strcmp(RegressFit,'Robust')
	[B,stats,ychap] = myRobustFit(x,y);
	%.........................
	xmin=min(x);
	xmax=max(x);
	dx=(xmax-xmin)/100;
	xxi=[xmin:dx:xmax];
% $$$ 	J=find(isnan(x)==0 & isnan(ychap)==0);
	J=find(isnan(x)==0 & isnan(ychap)==0 & x>0);
	yychap=ychap(J);
	xx=x(J);
	xchapI=xxi;
	ychapI=interp1(xx,yychap,xxi);
	%.........................
    end
end

%%%%%%
%PLOT:
%%%%%%
%................................
setXlim = [xzbitmap(1) xzbitmap(end)];
setYlim = [yzbitmap(1) yzbitmap(end)];
setXtick = xzbitmap(:);
setYtick = yzbitmap(:);
%................................
% $$$ setXtickLabel = xz;
% $$$ setYtickLabel = yz;
%................................
for i=1:length(xz)
    setXtickLabel{i} = xz(i);
end
%................................
for i=1:length(yz)
    setYtickLabel{i} = yz(i);
end
%................................
if strcmp(keyXdataNegative,'yes')
    %DATA:
    log10x = log10x*(-1);
    if strcmp(RegressLine,'yes')
	xchapI = xchapI*(-1);
    end
    %TICKS:
    %............
    xz = xz*(-1);
    xzbitmap = xzbitmap*(-1);
    %............
    setXlim = [xzbitmap(end) xzbitmap(1)];
    setXtick = flipud(xzbitmap(:));
    %............
% $$$     setXtickLabel = flipud(xz(:));
    %............
    for i=1:length(xz)
	xxz=flipud(xz(:));
	setXtickLabel{i} = xxz(i);
    end
    %............
end
%??????????????????????????????????????????
% $$$ log10x
% $$$ setXlim
% $$$ setXtick
% $$$ setXtickLabel
% $$$ pause
%??????????????????????????????????????????

%................................
%For modelNPHC resilience:
% $$$ %%setXtick(5)=[];
setXtickLabel{5}='';
%................................
% $$$ hp1=plot(log10x,log10y,'.'); %Plot all dots equal.
%................................
hp1=myPlot(log10x,log10y,'.'); %For different symbols (only works if x and y are arrays, otherwise all dots will be equal).
%................................
set(gca,'Xlim',setXlim,'Xtick',setXtick,'Xticklabel',setXtickLabel)
set(gca,'Ylim',setYlim,'Ytick',setYtick,'Yticklabel',setYtickLabel)
grid on
%................................
if strcmp(RegressLine,'yes')
    hold on
    hp2=plot(xchapI,ychapI,'k-'); %Add regression line.
    hold off
    set(gca,'Xlim',setXlim,'Xtick',setXtick,'Xticklabel',setXtickLabel)
    set(gca,'Ylim',setYlim,'Ytick',setYtick,'Yticklabel',setYtickLabel)
    grid on
else
    hp2=nan;
% $$$     set(gca,'Xlim',[xzbitmap(1) xzbitmap(end)],'Xticklabel','')
% $$$     set(gca,'Ylim',[yzbitmap(1) yzbitmap(end)],'Yticklabel','')
    %%grid on
end
%................................
