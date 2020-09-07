function []=myglobalmap(ANAME,A,Aminmax,cmap,pcnt,varargin)

%************************************
%Programa MYGLOBALMAPS.m:
%Use: 
%
% [] = myglobalmap('ANAME',A,Amimmax,cmap,pcnt)
%
% pcnt = define los limites sup-inf: 
%   limsup=maxA+(pcnt*maxA);
%   liminf=minA+(pcnt*minA);
% (usar por ej. 0.02, 0.035, etc. Depende del rango de valores del colormap.
%************************************

%%%%%%%%%%%%%%%%%%%
%OBTAIN ARRAY SIZE:
%%%%%%%%%%%%%%%%%%%
[msize,nsize]=size(A);

%%%%%%%%%%%%%%%%%
%OBTENGO LA LAND:
%%%%%%%%%%%%%%%%%
%===================================
[GLOBE,Land]=myglobeland(msize,nsize); %Land.
%===================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE IMAGESC LIMIT VALUES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================
if isempty(Aminmax)==0 %if NOT empty.
    %...............................
    minA=Aminmax(1);
    maxA=Aminmax(2);
    %...............................
else
    %...............................
    minA=min(A(:));
    maxA=max(A(:));
    %...............................
end
%...................................
Imin = find(A < minA);
Imax = find(A > maxA);
%...................................
A(Imin) = minA;
A(Imax) = maxA;
%...................................
Amax = max(A(:));
Amin = min(A(:));
%......................................
AmaxSign = sign(Amax);
AminSign = sign(Amin);
%......................................
if AminSign == 0
    AminSign = -1;
end
%......................................
Alimsup = (abs(Amax) + pcnt*abs(Amax))*AmaxSign
Aliminf = (abs(Amin) + pcnt*abs(Amax))*AminSign
%......................................
%===================================
%.......................................
if msize==180 %el mapa va de 90N a 90S.
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
elseif msize==162 %el mapa va de 80N a 80S.
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
myYticks = lattick3;
myXticks = longtick3;
myYticksLabel = gradoslat3;
myXticksLabel = gradoslong3;
%......................................
%===================================
%......................................
mySizeFontTicks = [10]; %for PreySwitching (4 panels).
myFontSizeTitle = [12];
%......................................
%===================================
%...................................
keyColorbarPlace='horiz';
% $$$ keyColorbarPlace='vertic';
%...................................
keyImageType='imagesc';
% $$$ keyImageType='contourf';
%...................................
% $$$ keyImageScale='Log2';
% $$$ keyImageScale='Log10';
keyImageScale='Linear';
%...................................
%===================================

%%%%%%%%%%%%%%%%%%%%
%ADD BLACK AND GRAY:
%%%%%%%%%%%%%%%%%%%%
%===================================
%...................................
MAP = cmap;
MAP = [[0 0 0];cmap;[0.5 0.5 0.5]];
%...................................
%===================================
%......................................
nbins = length(cmap);
% $$$ nbins = length(MAP);
%......................................
Adel = (Amax - Amin)/(nbins);
%......................................
hcbarlim = [Amin,Amax];
hcbartick = [Amin:Adel:Amax];
hcbarticklabel = hcbartick;
%......................................
% $$$ hcbarlim = [(Amin-Adel),(Amax+Adel)];
% $$$ hcbartick = [(Amin-Adel):Adel:(Amax+Adel)];
% $$$ hcbarticklabel = hcbartick;
%......................................
% $$$ hcbarlim = [0.1,0.9];
% $$$ hcbartick = [0.1:Adel:0.9];
% $$$ hcbarticklabel = hcbartick;
%......................................
%===================================
%......................................
Aimg = A;
Aimg(Land) = Alimsup;
%......................................
if strcmp(keyImageType,'imagesc')
    himg = imagesc(Aimg,[Aliminf Alimsup]);
    hcmap = colormap(MAP);
    hcbar = mycolorbar(keyColorbarPlace); %USAR ESTA.
    htitle=title(ANAME);
    caxis([Amin Amax])
    set(gca, 'CLim', [Amin, Amax])
    set(gca,'YTick',myYticks,'YTickLabel',myYticksLabel,'XTick',myXticks,'XTickLabel',myXticksLabel,'FontSize',mySizeFontTicks)
    set(hcbar ,'FontSize',mySizeFontTicks);
    set(htitle,'FontSize',myFontSizeTitle);
    if strcmp(keyImageScale,'Log10') | strcmp(keyImageScale,'Log2')
	if strcmp(keyColorbarPlace,'horiz')
	    set(hcbar,'Xlim',[zbitmap(1) zbitmap(end)],'Xtick',zbitmap(Itick),'Xticklabel',zdata(Itick))
	elseif strcmp(keyColorbarPlace,'vertic')
	    set(hcbar,'Ylim',[zbitmap(1) zbitmap(end)],'Ytick',zbitmap(Itick),'Yticklabel',zdata(Itick))
	end
    else
	if strcmp(keyColorbarPlace,'horiz')
	    set(hcbar,'Xlim',hcbarlim,'XTick',hcbartick,'XTickLabel',hcbarticklabel)
	elseif strcmp(keyColorbarPlace,'vertic')
	    set(hcbar,'Ylim',hcbarlim,'YTick',hcbartick,'YTickLabel',hcbarticklabel)
	end
    end
else
end
%===================================
