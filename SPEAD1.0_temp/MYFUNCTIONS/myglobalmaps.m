function [varargout]=myglobalmaps(VAR,VARNAME,VARminmax,cmap,pcnt,mode,varargin)
%************************************
%Programa MYGLOBALMAPS.m:
%Este programa grafica 12 subplots con monthly-means values.
%Uso: 
%
%[]=myglobalmaps(VAR,'VARNAME',cmap,pcnt,'mode',fignum,'scalemode','colorbarASK',cbarlims)
%
%mode='Monthly','Seasonal','Annual'.
%pcnt=define los limites sup-inf: 
%   limsup=maxVAR+(pcnt*maxVAR);
%   liminf=minVAR+(pcnt*minVAR);
%(usar por ej. 0.02, 0.035, etc. Depende del rango de valores del colormap.
%scalemode: 'Linear', 'Log10'.
%colorbarASK='yes', 'no'.
%cbarlims = vector de 3 componentes [cbartickinf,cbarticksup,cbardeltatick]
%************************************
%===================================
%...................................
VARNAME %tiene que ser un 'string'
%...................................
meses=['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
%...................................
%===================================
[m,n,p]=size(VAR)
% $$$ pcnt=0.020;
% $$$ %pcnt=0.025;
% $$$ %pcnt=0.035;
% $$$ % $$$ pcnt=0.05;
%...................................
%===================================
%OBTENGO LA LAND:
[GLOBE,Land]=myglobeland(m,n); %Land.
%...................................
%===================================
if isempty(VARminmax)==0 %if NOT empty.
    %...............................
    minVAR=VARminmax(1);
    maxVAR=VARminmax(2);
    %...............................
else
    %...............................
    minVAR=min(VAR(:));
    maxVAR=max(VAR(:));
    %...............................
end
%...................................
Imin = find(VAR < minVAR);
Imax = find(VAR > maxVAR);
%...................................
VAR(Imin) = minVAR;
VAR(Imax) = maxVAR;
%...................................
%===================================
for k=1:p
    %...............................
    VARk=VAR(:,:,k);
    %...............................
    VARk(end,end  )=maxVAR;
    VARk(end,end-1)=minVAR;
    %...............................
    VARout(:,:,k)=VARk;
    %...............................
end
VAR=VARout;
%===================================
%...................................
scalemode='Linear';
%...................................
% $$$ colorbartype='vertic';
colorbartype='horiz';
%...................................
colorbarASK='yes';
% $$$ colorbarASK='not';
%...................................
%===================================
%...................................
% $$$ mySizeFontTicks = [10]; %for PreySwitching (4 panels).
% $$$ myFontSizeTitle = [12];
%...................................
mySizeFontTicks = [8]; %for PreySwitching (16 panels).
myFontSizeTitle = [10];
%...................................
% $$$ mySizeFontTicks = [6]; %for PreySwitching (16 panels).
% $$$ myFontSizeTitle = [8];
%...................................
%===================================
% $$$ fignum=1;
nvarargin=length(varargin);
if nvarargin==1
    fignum=varargin{1};
elseif nvarargin==2
    fignum=varargin{1};
    scalemode=varargin{2};
    if strcmp(scalemode,'Log10') | strcmp(scalemode,'Log2')
	%.............................................
	Xdata = VAR;
	Xlims = VARminmax;
	%.............................................
	keyLogBase=scalemode;
	%.............................................
	[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup,zlims,Itick]=mybitmapscale(Xdata,Xlims,keyLogBase); %default.
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.0001);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.001);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.01);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.05);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.1);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.1,200); %for DMS.
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.05,200);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.05,1);
	%.............................................
	VAR=Xstar;
    end
elseif nvarargin==3
    fignum=varargin{1};
    scalemode=varargin{2};
    colorbarASK=varargin{3};
    if strcmp(scalemode,'Log10') | strcmp(scalemode,'Log2')
	%.............................................
	Xdata = VAR;
	Xlims = VARminmax;
	%.............................................
	keyLogBase=scalemode;
	%.............................................
	[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup,zlims,Itick]=mybitmapscale(Xdata,Xlims,keyLogBase); %default.
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.0001);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.001);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.01);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.05);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.1);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.1,200); %for DMS.
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.05,200);
	%.............................................
	%[Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,keyLogBase,0.05,1);
	%.............................................
	VAR=Xstar;
    end
elseif nvarargin==4
    fignum=varargin{1};
    scalemode=varargin{2};
    colorbarASK=varargin{3};
    cbarticklims=varargin{4};
    cbartickinf=cbarticklims(1);
    cbarticksup=cbarticklims(2);
    cbardtick=cbarticklims(3);
end
imagetype='imagesc';
%imagetype='contourf';
scalemode
colorbarASK
%===================================
%.......................................
if m==180 %el mapa va de 90N a 90S.
    %...............................................................
    gradoslat1=[+90,+80,+60,+40,+20,0,-20,-40,-60,-80,-90];
    gradoslong1=[-180,-160,-140,-120,-100,-80,-60,-40,-20,0,+20,+40,+60,+80,+100,+120,+140,+160,+180];
    lattick1=[1,[10:20:170],180];
    longtick1=[1,[20:20:360]];
    %...............................................................
    gradoslat2=[+80,+60,+40,+20,0,-20,-40,-60,-80];
    gradoslong2=[-160,-120,-80,-40,0,+40,+80,+120,+160];
    lattick2=[10:20:170];
    longtick2=[20:40:340];
    %...............................................................
    gradoslat3  = [+ 90 ,+60,+30, 00, -30, -60,-90];
    gradoslong3 = [-180,-120,-60, 00, +60,+120,+180];
    lattick3  = [01,[30:30:160],180];
    longtick3 = [01,[60:60:300],360];
    %...............................................................
    gradoslat4  = [+60,+30, 00, -30, -60];
    gradoslong4 = [-120,-60, 00, +60,+120];
    lattick4  = [30:30:160];
    longtick4 = [60:60:300];
    %...............................................................
elseif m==162 %el mapa va de 80N a 80S.
    gradoslat2=[+80,+60,+40,+20,0,-20,-40,-60,-80];
    gradoslong2=[-160,-120,-80,-40,0,+40,+80,+120,+160];
    lattick2=[1:20:170];
    longtick2=[20:40:340];
end
%......................................
%===================================
%......................................
% $$$ gradoslat = gradoslat2;
% $$$ gradoslong = gradoslong2;
% $$$ lattick = lattick2;
% $$$ longtick = longtick2;
%......................................
% $$$ myYticks = lattick3;
% $$$ myXticks = longtick3;
% $$$ myYticksLabel = gradoslat3;
% $$$ myXticksLabel = gradoslong3;
%......................................
myYticks = lattick4;
myXticks = longtick4;
myYticksLabel = gradoslat4;
myXticksLabel = gradoslong4;
%......................................
%===================================
%......................................
maxVAR=max(VAR(:));
minVAR=min(VAR(:));
%......................................
limsup=maxVAR+(pcnt*maxVAR);
%......................................
if minVAR>=0
    liminf=minVAR-(pcnt*maxVAR);
elseif minVAR<0
    liminf=minVAR-(pcnt*abs(minVAR));
    liminf=-max([abs(limsup),abs(liminf)]);
end
%......................................
% $$$ limsup = ceil(limsup*100)/100;
% $$$ liminf = floor(liminf*1000)/1000;
%......................................
minmaxVAR=[minVAR,maxVAR]
liminfsup=[liminf,limsup]
%......................................
hcbarlim=[minVAR,maxVAR];
if nvarargin<4
    %......................................
% $$$     nbins=10;
% $$$     nbins=5;
    nbins=4;
    %......................................
% $$$     cbardtick=(limsup+abs(liminf))/nbins;
% $$$     hcbartick=[liminf:cbardtick:limsup];
    %......................................
% $$$     cbardtick = limsup/nbins;
% $$$     hcbartick = [0:cbardtick:limsup];
    %......................................
    cbardtick = maxVAR/nbins;
    hcbartick = [0:cbardtick:maxVAR];
    %......................................
    hcbarticklabel=hcbartick;
elseif nvarargin==4
    hcbartick=[cbartickinf:cbardtick:cbarticksup];
    hcbarticklabel=hcbartick;
end
hcbarlim
hcbartick
hcbarticklabel
%......................................
%%%%%%%%%%%%%%%%%%%%
%ADD BLACK AND GRAY:
%%%%%%%%%%%%%%%%%%%%
whos cmap
% $$$ MAP = cmap; %Use this one alone para Fig0 de PNASpaper (i.e. Do NOT use [0 0 0; MAP; 0.5 0.5 0.5]).
MAP = [[0 0 0];cmap;[0.5 0.5 0.5]];
%......................................
% $$$ if strcmp(colorbarASK,'not')
% $$$     colorbartype='delete';
% $$$ end
%......................................
if strcmp(class(fignum),'char') %si fignum='junk'
else
    figure(fignum)
end
%......................................
ii=3;jj=4;
if strcmp(mode,'Monthly')
    if p==12
	ii=3;jj=4;
    elseif p==6
	ii=3;jj=2;
    elseif p==4
	ii=2;jj=2;
    elseif p==2
	ii=1;jj=2;
    elseif p==72
	ii=8;jj=9;
    end
    for k=1:p
	subplot(ii,jj,k)
	VARk=VAR(:,:,k);
	VARk(Land)=limsup;
	VARk(1,1)=liminf;
	VARk(1,2)=limsup;
	%..................................
	imagesc(VARk);
	colormap(MAP)
	hcbar=mycolorbar(colorbartype);
	%..................................
	%Lo repito para que se ponga bien la "colorbar" (hay un bug cuando uso
	%mis propios "colormaps").
	imagesc(VARk);
	colormap(MAP)
	hcbar=mycolorbar(colorbartype);
	%..................................
	if p==12
	    mesk=meses(k,:);
	    htitle=title([VARNAME,' ',mesk]);
	    set(htitle,'FontSize',mySizeFontTicks)
	end
	set(gca,'YTick',myYticks,'YTickLabel',myYticksLabel,'XTick',myXticks,'XTickLabel',myXticksLabel,'FontSize',mySizeFontTicks)
	set(hcbar,'FontSize',mySizeFontTicks);
	if strcmp(scalemode,'Log10') | strcmp(scalemode,'Log2')
	    if strcmp(colorbartype,'horiz')
		set(hcbar,'Xlim',[zbitmap(1) zbitmap(end)],'Xtick',zbitmap(:),'Xticklabel',zdata)
	    elseif strcmp(colorbartype,'vertic')
		set(hcbar,'Ylim',[zbitmap(1) zbitmap(end)],'Ytick',zbitmap(:),'Yticklabel',zdata)
	    end
	end
	grid on
    end
elseif strcmp(mode,'Seasonal')
    Iwinter=[12,1,2];
    Ispring=[3,4,5];
    Isummer=[6,7,8];
    Iautom=[9,10,11];

    %......................................
% $$$     VARwinter=snanmean(VAR(:,:,Iwinter),3);
% $$$     VARspring=snanmean(VAR(:,:,Ispring),3);
% $$$     VARsummer=snanmean(VAR(:,:,Isummer),3);
% $$$     VARautom=snanmean(VAR(:,:,Iautom),3);
    %......................................
    VARwinter=nanmean(VAR(:,:,Iwinter),3);
    VARspring=nanmean(VAR(:,:,Ispring),3);
    VARsummer=nanmean(VAR(:,:,Isummer),3);
    VARautom =nanmean(VAR(:,:,Iautom ),3);
    %......................................

% $$$     if strcmp(VARNAME,'ETA')==0 %si uso ETA no pongo Land mask (para el resto de variables si).
    VARwinter(Land)=limsup;
    VARspring(Land)=limsup;
    VARsummer(Land)=limsup;
    VARautom(Land)=limsup;
% $$$     end    

% $$$     VARwinter(1,1)=liminf;VARwinter(1,2)=limsup;
% $$$     VARspring(1,1)=liminf;VARspring(1,2)=limsup;
% $$$     VARsummer(1,1)=liminf;VARsummer(1,2)=limsup;
% $$$     VARautom(1,1)=liminf;VARautom(1,2)=limsup;
    
    subplot(2,2,1)
    imagesc(VARwinter,[liminf limsup])
    htitle=title([VARNAME,' ','DJF']); % (austral-summer)
    hcbar=mycolorbar(colorbartype);
    colormap(MAP)
    %..................................
    %Lo repito para que se ponga bien la "colorbar" (hay un bug cuando uso
    %mis propios "colormaps").
    hcbar=mycolorbar(colorbartype);
    colormap(MAP)
    %..................................
    set(gca,'YTick',myYticks,'YTickLabel',myYticksLabel,'XTick',myXticks,'XTickLabel',myXticksLabel,'FontSize',mySizeFontTicks)
    set(hcbar,'FontSize',mySizeFontTicks);
    set(htitle,'FontSize',myFontSizeTitle)
    if strcmp(scalemode,'Log10') | strcmp(scalemode,'Log2')
	if strcmp(colorbartype,'horiz')
	    set(hcbar,'Xlim',[zbitmap(1) zbitmap(end)],'Xtick',zbitmap(:),'Xticklabel',zdata)
	elseif strcmp(colorbartype,'vertic')
	    set(hcbar,'Ylim',[zbitmap(1) zbitmap(end)],'Ytick',zbitmap(:),'Yticklabel',zdata)
	end
    end
    grid on

    subplot(2,2,2)
    imagesc(VARspring,[liminf limsup])
    htitle=title([VARNAME,' ','MAM']);% (austral-autom)
    hcbar=mycolorbar(colorbartype);
    colormap(MAP)
    set(gca,'YTick',myYticks,'YTickLabel',myYticksLabel,'XTick',myXticks,'XTickLabel',myXticksLabel,'FontSize',mySizeFontTicks)
    set(hcbar,'FontSize',mySizeFontTicks);
    set(htitle,'FontSize',myFontSizeTitle)
    if strcmp(scalemode,'Log10') | strcmp(scalemode,'Log2')
	if strcmp(colorbartype,'horiz')
	    set(hcbar,'Xlim',[zbitmap(1) zbitmap(end)],'Xtick',zbitmap(:),'Xticklabel',zdata)
	elseif strcmp(colorbartype,'vertic')
	    set(hcbar,'Ylim',[zbitmap(1) zbitmap(end)],'Ytick',zbitmap(:),'Yticklabel',zdata)
	end
    end
    grid on

    subplot(2,2,3)
    imagesc(VARsummer,[liminf limsup])
    htitle=title([VARNAME,' ','JJA']);% (austral-winter)
    hcbar=mycolorbar(colorbartype);
    colormap(MAP)
    set(gca,'YTick',myYticks,'YTickLabel',myYticksLabel,'XTick',myXticks,'XTickLabel',myXticksLabel,'FontSize',mySizeFontTicks)
    set(hcbar,'FontSize',mySizeFontTicks);
    set(htitle,'FontSize',myFontSizeTitle)
    if strcmp(scalemode,'Log10') | strcmp(scalemode,'Log2')
	if strcmp(colorbartype,'horiz')
	    set(hcbar,'Xlim',[zbitmap(1) zbitmap(end)],'Xtick',zbitmap(:),'Xticklabel',zdata)
	elseif strcmp(colorbartype,'vertic')
	    set(hcbar,'Ylim',[zbitmap(1) zbitmap(end)],'Ytick',zbitmap(:),'Yticklabel',zdata)
	end
    end
    grid on
 
    subplot(2,2,4,[liminf limsup])
    imagesc(VARautom)
    htitle=title([VARNAME,' ','SON']);% (austral-spring)
    hcbar=mycolorbar(colorbartype);
    colormap(MAP)
    set(gca,'YTick',myYticks,'YTickLabel',myYticksLabel,'XTick',myXticks,'XTickLabel',myXticksLabel,'FontSize',mySizeFontTicks)
    set(hcbar,'FontSize',mySizeFontTicks);
    set(htitle,'FontSize',myFontSizeTitle)
    if strcmp(scalemode,'Log10') | strcmp(scalemode,'Log2')
	if strcmp(colorbartype,'horiz')
	    set(hcbar,'Xlim',[zbitmap(1) zbitmap(end)],'Xtick',zbitmap(:),'Xticklabel',zdata)
	elseif strcmp(colorbartype,'vertic')
	    set(hcbar,'Ylim',[zbitmap(1) zbitmap(end)],'Ytick',zbitmap(:),'Yticklabel',zdata)
	end
    end
    grid on
    
% $$$ elseif strcmp(mode,'Biannual')
% $$$     Ispringsummer=[4:9];
% $$$     Iautumwinter=[10,11,12,1,2,3];

elseif strcmp(mode,'Annual')
    if p>1 
	disp('Error: VAR tiene que ser matriz 2D!!')
    end
    %..................................
    VVAR=VAR;
    %..................................
    if p==1 %ONLY FOR 2D ARRAYS.
	VVAR(Land)=limsup; %(Uncomment for RGB plot of GlobalPDRpaper)
    end
    %..................................
    VVAR(1,1)=liminf;
    VVAR(1,2)=limsup;
    %..................................
% $$$     VVAR=VAR; %para Fig0 de PNASpaper.
% $$$     Inan=find(isnan(VAR));
% $$$     VVAR(Inan)=-1-0.08; %NEGRO(RRsig = pixels con R significativa (el resto es -1.01) 
% $$$     VVAR(Land)=-1-0.04; %GRIS
    %..................................
    VVARminmax=[min(VVAR(:)),max(VVAR(:))]
    if strcmp(imagetype,'imagesc')
	%..................................
% $$$ 	himg=imagesc(VVAR);
% $$$ 	%..................................
	himg=imagesc(VVAR,[liminf limsup])
% $$$ 	%..................................
% $$$ 	himg=imagesc(VVAR,[minVAR maxVAR])
	%..................................
% $$$ 	set(gca,'Color',[0 0 0]);
% $$$ 	set(himg,'AlphaData',0.5);
	%..................................
    elseif strcmp(imagetype,'contourf')
% $$$     Inan=find(isnan(VVAR)==1);
% $$$     VVAR(Inan)=liminf;
% $$$ 	[C,h]=contourf(flipud(VVAR),[liminf,[-1:0.2:+1],limsup]);
	[C,h]=contourf(flipud(VVAR),[-1:0.2:+1]);
    end
    %..................................
    hcmap = colormap(MAP);
    %..................................
    if strcmp(colorbarASK,'yes')
% $$$     hcbar = colorbar(colorbartype); %For RGB map (globalPDRpaper).
    hcbar = mycolorbar(colorbartype); %USAR ESTA.
    elseif strcmp(colorbarASK,'not')
	hcbar = [];
    end
    %..................................
    %ALLOW DIFFERENT COLORMAPS FOR THE SUBPLOTS:
% $$$     pause(1)
% $$$     freezeColors
% $$$     cbfreeze(hcbar)
    %..................................

% $$$     %Lo repito para que se ponga bien la "colorbar" (hay un bug cuando uso mis propios "colormaps").
% $$$     if strcmp(imagetype,'imagesc')
% $$$ 	himg=imagesc(VVAR);
% $$$     elseif strcmp(imagetype,'contourf')
% $$$ % $$$ 	[C,h]=contourf(flipud(VVAR),[liminf,[-1:0.2:+1],limsup]);
% $$$ 	[C,h]=contourf(flipud(VVAR),[-1:0.2:+1]);
% $$$ 	clabel(C,h,[0.4],'LabelSpacing',100)
% $$$     end
    %..................................
% $$$     hcmap = colormap(MAP);
% $$$     hcbar = mycolorbar(colorbartype); %USAR ESTA.
    %..................................
% $$$     hcbarPositionOld  = get(hcbar,'Position');
% $$$     hcbarPositionNew(1) = hcbarPositionOld(1) + 0.00*hcbarPositionOld(1);
% $$$     hcbarPositionNew(2) = hcbarPositionOld(2) + 0.00*hcbarPositionOld(2);
% $$$     hcbarPositionNew(3) = hcbarPositionOld(3) + 0.00*hcbarPositionOld(3);
% $$$     hcbarPositionNew(4) = hcbarPositionOld(4) + 0.00*hcbarPositionOld(4);
% $$$     set(hcbar,'Position',hcbarPositionNew);
    %..................................
    %============================================================================
    %----------------------------------------------------------------------------
    %NOTE: pos = [left, bottom, width, height]
    %----------------------------------------------------------------------------
    %For PreySwitchingPaper 16 plots:
% $$$     axesPosition = get(gca,'Position');
% $$$     hcbar = colorbar('location','South','XAxisLocation','bottom')
% $$$     hcbarPositionOld  = get(hcbar,'Position');
% $$$     hcbarPositionNew(1) = axesPosition(1);
% $$$     hcbarPositionNew(2) = hcbarPositionOld(2) - 2.00*hcbarPositionOld(4);
% $$$     hcbarPositionNew(3) = axesPosition(3);
% $$$     hcbarPositionNew(4) = hcbarPositionOld(4) - 0.25*hcbarPositionOld(4);
% $$$     set(hcbar,'Position',hcbarPositionNew);
    %..................................
% $$$     pause
% $$$     axes1 = gca;
% $$$     axes1position = get(axes1, 'Position');
% $$$     axes2 = axes('Position', axes1position,'Visible', 'off');
% $$$     hcbar = mycolorbar(colorbartype,'peer',axes2);
    %..................................
    %============================================================================
    htitle=title(VARNAME);
    %...........................
% $$$     set(gca,'YTick',myYticks,'YTickLabel',[],'XTick',myXticks,'XTickLabel',[],'FontSize',mySizeFontTicks)
% $$$     set(gca,'YTick',myYticks,'YTickLabel',[],'XTick',myXticks,'XTickLabel',myXticksLabel,'FontSize',mySizeFontTicks)
    %...........................
    set(gca,'YTick',myYticks,'YTickLabel',myYticksLabel,'XTick',myXticks,'XTickLabel',myXticksLabel,'FontSize',mySizeFontTicks)
% $$$     set(gca,'YTick',myYticks,'YTickLabel',myYticksLabel,'XTick',myXticks,'XTickLabel',[],'FontSize',mySizeFontTicks)
    %...........................
% $$$     if strcmp(colorbartype,'horiz')
% $$$ 	if minVAR==-1 & maxVAR==+1 %para "corrGLOBAL.m".
% $$$ 	    set(hcbar,'Xlim',[-1 +1],'XTick',[-1:0.2:+1],'XTickLabel',[-1:0.2:+1])
% $$$ 	elseif minVAR>-0.9 & maxVAR<+0.9 %para "dmsalgGW.m".
% $$$ 	    set(hcbar,'Xlim',[-0.8 +0.8],'XTick',[-0.8:0.2:+0.8],'XTickLabel',[-0.8:0.2:+0.8])
% $$$ 	elseif minVAR>-35 & maxVAR>=+15 & maxVAR<+35 %para "dmsalgGW.m".
% $$$ 	    set(hcbar,'Xlim',[-30 +30],'XTick',[-30:5:+30],'XTickLabel',[-30:5:+30])
% $$$ 	elseif minVAR>=0
% $$$ 	    set(hcbar,'Xlim',[0 maxVAR],'XTick',[0:1:maxVAR],'XTickLabel',[0:1:maxVAR])
% $$$ 	end
% $$$     elseif strcmp(colorbartype,'vertic')
% $$$ 	if minVAR==-1 & maxVAR==+1 %para "corrGLOBAL.m".
% $$$ 	    set(hcbar,'Ylim',[-1 +1],'YTick',[-1:0.2:+1],'YTickLabel',[-1:0.2:+1])
% $$$ 	elseif minVAR>-0.9 & maxVAR<+0.9 %para "dmsalgGW.m".
% $$$ 	    set(hcbar,'Ylim',[-0.8 +0.8],'YTick',[-0.8:0.2:+0.8],'YTickLabel',[-0.8:0.2:+0.8])
% $$$ 	elseif minVAR>-35 & maxVAR>=+15 & maxVAR<+35 %para "dmsalgGW.m".
% $$$ 	    set(hcbar,'Ylim',[-30 +30],'YTick',[-30:5:+30],'YTickLabel',[-30:5:+30])
% $$$ 	elseif minVAR>=0
% $$$ 	    set(hcbar,'Ylim',[0 maxVAR],'YTick',[0:1:maxVAR],'YTickLabel',[0:1:maxVAR])
% $$$ 	end
% $$$     end
    %...........................
    if strcmp(colorbarASK,'yes')
    if strcmp(scalemode,'Log10') | strcmp(scalemode,'Log2')
% $$$ 	mycmapZlims = [zbitmap(1) zbitmap(end)];
	mycmapZlims = [zlims(1) zlims(end)];
	if strcmp(colorbartype,'horiz')
	    set(hcbar,'Xlim',mycmapZlims,'Xtick',zbitmap(Itick),'Xticklabel',zdata(Itick))
	elseif strcmp(colorbartype,'vertic')
	    set(hcbar,'Ylim',mycmapZlims,'Ytick',zbitmap(Itick),'Yticklabel',zdata(Itick))
	end
    else
	if strcmp(colorbartype,'horiz')
	    set(hcbar,'Xlim',hcbarlim,'XTick',hcbartick,'XTickLabel',hcbarticklabel)
	elseif strcmp(colorbartype,'vertic')
	    set(hcbar,'Ylim',hcbarlim,'YTick',hcbartick,'YTickLabel',hcbarticklabel)
	end
    end
    end
    %...........................
    grid on
    %....................................
    if strcmp(colorbarASK,'yes')
    set(hcbar ,'FontSize',mySizeFontTicks);
    end
    set(htitle,'FontSize',myFontSizeTitle);
% $$$     set(htitle,'Visible','off') %For PreySwitching paper 16 panels.
    %....................................
% $$$     if strcmp(colorbarASK,'not')
% $$$ 	if strcmp(colorbartype,'horiz')
% $$$ 	    delete(colorbar('horiz'))
% $$$ 	elseif strcmp(colorbartype,'vertic')
% $$$ 	    delete(colorbar('vertic'))
% $$$ 	end
% $$$     end
    %....................................
    if strcmp(colorbarASK,'not')
	set(hcbar,'Visible','off')
% $$$     colorbar('off')
% $$$     colorbar('hide')
% $$$     colorbar('delete')
    end
    %....................................
    set(gca,'YAxisLocation','left') %Yticks on the left side.
% $$$     set(gca,'YAxisLocation','right') %Yticks on the right side.
    %....................................
end

%%%%%%%%
%OUTPUT:
%%%%%%%%
varargout{1} = himg;
varargout{2} = hcbar;
