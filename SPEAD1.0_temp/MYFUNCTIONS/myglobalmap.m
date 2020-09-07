function [hcbar] = myglobalmap(ANAME,A,Aminmax,cmap,pcnt,keyImageScale)

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
%OBTAIN A GIVEN PERCENTILE OF THE DATA THAT ARE ABOVE THE MEAN:
%...................................
B = A(:);
%...................................
B = B(B > 0); %Use only positive values (no zero or nan)
%...................................
B(B < mean(B)) = []; %Remove points below mean average.
%...................................
A50mean = B; %Only data points above mean value.
%...................................
%%Pcntile = 80; %Percentile 80%
Pcntile = 99.9; %Percentile 80%
%...................................
A80pcntile = prctile(B,Pcntile);
%...................................
%===================================
if isempty(Aminmax)==0 %if NOT empty.
    %...............................
    minA = Aminmax(1);
    maxA = Aminmax(2);
    %...............................
else %if Aminmax is empty.
    %...............................
    minA = min(A(:));
    maxA = max(A(:));
    %...............................
% $$$     minA = min(A(:));
% $$$     maxA = A80pcntile; %Percentile 80% of the data points above mean value.
    %...............................
    Aminmax = [minA,maxA];
    %...............................
end
%...................................
Imin = find(A < minA);
Imax = find(A > maxA);
%...................................
A(Imin) = minA;
A(Imax) = maxA;
%...................................
% $$$ Amax = max(A(:)); %Optional (comment-out otherwise)
% $$$ Amin = min(A(:));
%......................................
Amax = maxA; %Optional (comment-out otherwise)
Amin = minA;
%......................................
AmaxSign = sign(Amax);
AminSign = sign(Amin);
%......................................
% $$$ if AminSign == 0
% $$$     AminSign = -1;
% $$$ end
%......................................
AminSign = -1;
%......................................
% $$$ Alimsup = (abs(Amax) + pcnt*abs(Amax))*AmaxSign;
% $$$ Aliminf = (abs(Amin) + pcnt*abs(Amax))*AminSign;
%......................................
Alimsup = abs(Amax) + (pcnt*abs(Amax))*AmaxSign;
Aliminf = abs(Amin) + (pcnt*abs(Amax))*AminSign;
%......................................
% $$$ Aliminf = Amin;
%......................................
if strcmp(keyImageScale,'Logscale')


end
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
    gradoslat4  = [+60,+30, 00, -30, -60];
    gradoslong4 = [-120,-60, 00, +60,+120];
    lattick4  = [30:30:160];
    longtick4 = [60:60:300];
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
mySizeFontTicks = [10]; %for PreySwitching (4 panels).
myFontSizeTitle = [12];
%......................................
%===================================
%...................................
% $$$ % $$$ keyImageScale='Log2';
% $$$ % $$$ keyImageScale='Log10';
% $$$ keyImageScale='Linear';
%...................................
keyImageType='imagesc';
% $$$ keyImageType='contourf';
%...................................
%%keyColorbarPlace='vertic';
keyColorbarPlace='horiz';
%...................................
%===================================

%%%%%%%%%%
%LOGSCALE:
%%%%%%%%%%
if strcmp(keyImageScale,'Logscale')
    [Amap,zmap,zdat] = myBitmap256conversion(A,'Logscale');
    [Amapminmax] = myBitmap256conversion(Aminmax,'Logscale');
    A = Amap;
    Amin = Amapminmax(1);
    Amax = Amapminmax(2);
    Aliminf = Amin;
    Alimsup = Amax;
end

%%%%%%%%%%%%%%%%%%%%
%ADD BLACK AND GRAY:
%%%%%%%%%%%%%%%%%%%%
%===================================
%...................................
MAP = cmap;
%...................................
MAP = [cmap;[0.5 0.5 0.5]]; %Gray map for Land. 
%...................................
% $$$ MAP = [[0 0 0];cmap;[0.5 0.5 0.5]]; %Gray map for Land + Black map for no-data.
%...................................
%===================================
%......................................
nbins = 10;
% $$$ nbins = length(cmap);
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
himg = imagesc(Aimg);
set(gca,'CLim', [Aliminf, Alimsup]) %Place clim before the colorbar
caxis([Aliminf Alimsup])
hcmap = colormap(MAP);
%%hcbar = mycolorbar(keyColorbarPlace); 
%%hcbar = mycolorbarMatlab5p0(keyColorbarPlace); 
hcbar = mycolorbarMatlab6p5(keyColorbarPlace); 
htitle=title(ANAME);
set(gca,'YTick',myYticks,'YTickLabel',myYticksLabel)
set(gca,'XTick',myXticks,'XTickLabel',myXticksLabel)
set(gca,'FontSize',mySizeFontTicks)
set(hcbar ,'FontSize',mySizeFontTicks);
set(htitle,'FontSize',myFontSizeTitle);
if strcmp(keyImageScale,'Logscale')
    if strcmp(keyColorbarPlace,'horiz')
	set(hcbar,'Xtick',zmap,'Xticklabel',zdat)
    elseif strcmp(keyColorbarPlace,'vertic')
	set(hcbar,'Ytick',zmap,'Yticklabel',zdat)
    end
elseif strcmp(keyImageScale,'Linscale')
% $$$     if strcmp(keyColorbarPlace,'horiz')
% $$$ 	set(hcbar,'Xlim',hcbarlim,'XTick',hcbartick,'XTickLabel',hcbarticklabel)
% $$$     elseif strcmp(keyColorbarPlace,'vertic')
% $$$ 	set(hcbar,'Ylim',hcbarlim,'YTick',hcbartick,'YTickLabel',hcbarticklabel)
% $$$     end
end
grid on
%===================================
