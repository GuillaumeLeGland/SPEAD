function [hp1,hp2]=myRegressLog10line(xinput,yinput,varargin)

%***************************************************    
%Function: MYREGRESSLOG10LINE.m
%
%Use: [hp1,hp2]=myRegressLog10line(xinput,yinput)
%Use: [hp1,hp2]=myRegressLog10line(xinput,yinput,[],[],[],'yes')
%Use: [hp1,hp2]=myRegressLog10line(xinput,yinput,stdx,[],[],'yes')
%Use: [hp1,hp2]=myRegressLog10line(xinput,yinput,stdx,stdy,[],'yes')
%Use: [hp1,hp2]=myRegressLog10line(xinput,yinput,stdx,stdy,Wxy,'yes','Robust')
%
%Inputs:
%
%xinput      = Vector (or array) of X-variable.
%yinput      = Vector (or array) of Y-variable.
%xstd = 
%ystd =

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

nvarargin=length(varargin);
%...................
% $$$ RegressLine='yes';
% $$$ RegressLine='no';
%...................
RegressFit='Normal';
% $$$ RegressFit='Robust';
%...................
ErrorbarXYshadow='no';
% $$$ ErrorbarXYshadow='yes';
%...................
xmin=min(xinput(:));
xmax=max(xinput(:));
ymin=min(yinput(:));
ymax=max(yinput(:));
%...................
stdx=[];
stdy=[];
%...................
if nvarargin==1
    stdx=varargin{:,1};
elseif nvarargin==2
    stdx=varargin{:,1};
    stdy=varargin{:,2};
elseif nvarargin==3
    %...................
    stdx=varargin{:,1};
    stdy=varargin{:,2};
    %...................
    xyminmax=varargin{:,3}
    %...................
    xmin=xyminmax(1,1);
    xmax=xyminmax(1,2);
    ymin=xyminmax(2,1);
    ymax=xyminmax(2,2);
    %...................
elseif nvarargin==4
    %...................
    stdx=varargin{:,1};
    stdy=varargin{:,2};
    %...................
    xyminmax=varargin{:,3}
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
    RegressLine=varargin{:,4}
    %...................
elseif nvarargin==5
    %...................
    stdx=varargin{:,1};
    stdy=varargin{:,2};
    %...................
    xyminmax=varargin{:,3}
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
    RegressLine=varargin{:,4}
    RegressFit=varargin{:,5}
    %...................
end
%.........................
% $$$ if isempty(stdx)==0 %If not empty.
% $$$     ErrorbarXYshadow='yes';
% $$$ end
%.........................
ErrorbarXYshadow='no';
%.........................
xminmax=[xmin,xmax]
yminmax=[ymin,ymax]
%.........................
    
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
[Xstar,xz,xzbitmap,xzbitmapinf,xzbitmapsup]=mybitmapscale(xinput,xmin,xmax);
[Ystar,yz,yzbitmap,yzbitmapinf,yzbitmapsup]=mybitmapscale(yinput,ymin,ymax);

log10x=Xstar;
log10y=Ystar;

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
%============================================
if strcmp(ErrorbarXYshadow,'yes')
    %.................
    [stdXstar]=mybitmapscale(stdx);
    [stdYstar]=mybitmapscale(stdy);
    %.................
% $$$     figure(20)
% $$$     subplot(2,2,1)
% $$$     plot(stdx,stdXstar,'r*')
% $$$     subplot(2,2,2)
% $$$     plot(stdy,stdYstar,'r*')
    %.................
    fcx=stdx./xinput;
    fcy=stdy./yinput;
    log10stdx=stdXstar.*fcx;
    log10stdy=stdYstar.*fcy;
    %.................
    log10stdx
    log10stdy
    %.................
    [hpatch]=myErrorbarXYshadow(log10x,log10y,log10stdx/2,log10stdy/2);
    %.................
    pause
end
%============================================
hp1=plot(log10x,log10y,'.'); %Plot all dots.
%============================================
if strcmp(RegressLine,'yes')
    hold on
    hp2=plot(xchapI,ychapI,'k-'); %Add regression line.
    hold off
    %set(gca,'Xlim',[xzbitmap(1) xzbitmap(end)],'Xtick',xzbitmap(:),'Xticklabel',xz)
    %set(gca,'Ylim',[yzbitmap(1) yzbitmap(end)],'Ytick',yzbitmap(:),'Yticklabel',yz)
    grid on
else
    hp2=nan;
% $$$     set(gca,'Xlim',[xzbitmap(1) xzbitmap(end)],'Xticklabel','')
% $$$     set(gca,'Ylim',[yzbitmap(1) yzbitmap(end)],'Yticklabel','')
    %%grid on
end
%============================================
