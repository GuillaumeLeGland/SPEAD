function [ax1,ax2,hcbar001,hcbar002,varargout] = myimagesc2Dtimedepth(A,Aminmax,ColorbarPosition,jend,dj,zzend,dz,varargin)
    
%****************************************
%Use: myimagesc2Dtimedepth(A,ColorbarPosition,jend,dj,zzend,dz,varargin)
%ColorbarPosition: ['Horizontal'] ['Vertical']
%varargin{1}=showlabelY ('yes','not') %by default 'yes'
%varargin{2}=showlabelX ('yes','not') %by default 'yes'
%varargin{3}=scalefromzero ('yes','not') %by default 'yes'
%varargin{4}=sizefont %by default is 10.
%****************************************

%...................................
% $$$ AA=ver;
% $$$ keyNameLogiciel=AA(1).Name;
% $$$ if strcmp(keyNameLogiciel,'Octave')
% $$$     keyOctave='yes';
% $$$ else
% $$$     keyOctave='not';
% $$$ end
%...................................
mv=version;
keyNameLogiciel=mv(1:3);
if strcmp(keyNameLogiciel,'3.2')
    keyOctave='yes';
else
    keyOctave='not';
end
%...................................
if strcmp(ColorbarPosition,'Horizontal')
    if strcmp(keyOctave,'not')
	barpos='horiz';
    elseif strcmp(keyOctave,'yes')
	barpos='SouthOutside';
    end
elseif strcmp(ColorbarPosition,'Vertical')
    if strcmp(keyOctave,'not')
	barpos='vertic';
    elseif strcmp(keyOctave,'yes')
	barpos='EastOutside';
    end
end
%...................................
%%keyColorbar = 'not';
keyColorbar = 'yes'; 
%...................................
keyPcolor='not';
% $$$ keyPcolor='yes';
%...................................
% $$$ ndays=365;
ndays=360;
%...................................
% $$$ monthlim=[1,31,59,90,120,151,181,212,243,273,304,334,365];
%...................................
monthlim=[0:30:ndays];
%...................................
monthlim2=0.5*diff(monthlim)+monthlim(1:end-1);
xticklabels=['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
%...................................
nvarargin=length(varargin);
Amin=0;
if nvarargin==0
    showlabelY='yes';
    showlabelX='yes';
    scalefromzero='yes';
    sizefont=10;
elseif nvarargin==1
    showlabelY=varargin{1};
    showlabelX='yes';
    scalefromzero='yes';
    sizefont=10;
elseif nvarargin==2
    showlabelY=varargin{1};
    showlabelX=varargin{2};
    scalefromzero='yes';
    sizefont=10;
elseif nvarargin==3
    showlabelY=varargin{1};
    showlabelX=varargin{2};
    scalefromzero=varargin{3};
    sizefont=10;
elseif nvarargin==4
    showlabelY=varargin{1};
    showlabelX=varargin{2};
    scalefromzero=varargin{3};
    sizefont=varargin{4};
end
%...................................
%NOTE: Matlab does * NOT * accept independent SizeFonts for Xticks and Yticks.
SizeFontXticks = sizefont;
SizeFontYticks = sizefont;
SizeFontColorbar = sizefont;
%...................................
%.............
WHITE=ones(64,3)*0.8;
%.............
MAP=jet;
% $$$ MAP=gray;
%.............
% $$$ maxA=max(A(:));
% $$$ if maxA==-999
% $$$     MAP(1,:)=[1 1 1];
% $$$ end
%.............
% $$$ maxA=max(A(:));
% $$$ if maxA~=-999 %si no son todo -999.
% $$$     A(1,1)=0; %for the colorbar.
% $$$ end
%.............
if strcmp(scalefromzero,'yes')
    A(1,1)=0; %for the colorbar.
end
%...................................
%MAP(1,:)=[1 1 1]; %for the colormap.
%MAP(1,:)=[0.5 0.5 0.5]; %for the colormap.
pcnt=0.000;
%%pcnt=0.025; %To extend colorbar limits.
if isempty(Aminmax) == 0 %if *not* empty.
    Amin = Aminmax(1);
    Amax = Aminmax(2);
elseif isempty(Aminmax) == 1 %if empty.
    Amax = max(A(:)) + sqrt(eps);
    Amin = min(A(:));
end

limsup = Amax+(pcnt*Amax);
if Amin >= 0
    liminf = Amin-(pcnt*Amax);
elseif Amin < 0
    liminf = Amin-(pcnt*abs(Amin));
    liminf = -max([abs(limsup),abs(liminf)]);
end
A(1,2) = liminf;

%Pongo NaN donde hay -999:
I999=find(A==-999);
A(I999)=nan;
%.............
minA = min(A(:));
maxA = max(A(:));

Aminmax = [Amin,Amax]; %show this.
minmaxA = [minA,maxA];

%%%%%%%%%%
%1st AXES:
%%%%%%%%%%
%............
[mm,nn]=size(A);
%............
Atmp=ones(mm,nn)*nan;
%............
% $$$ if minA~=maxA
% $$$     himg=imagesc(A,[minA maxA]);
% $$$ else
% $$$     himg=imagesc(A);
% $$$ end
%............
if strcmp(keyPcolor,'not')
    %..............
    himg1=imagesc(A);
    %..............
    %%himg1=imagesc(Atmp);
    %..............
elseif strcmp(keyPcolor,'yes')
    himg1=pcolor(flipud(A));
    shading interp
end
%............
% $$$ colormap(MAP)
% $$$ colormap(WHITE)
%...............
%Ha de estar before de "colorbar".
ax1 = gca;
ax1position = get(ax1, 'Position');
%...............
if strcmp(keyOctave,'yes')
    hcbar001 = colorbar(barpos);
else
    hcbar001 = mycolorbar(barpos);
end
%...............
myXlims = [1 ndays];
myXtick = monthlim2;
myXtickLabel = xticklabels;
%...............
% $$$ myYlims   = [0 jend];
% $$$ myYtick = [0:dj:jend];
%...............
myYlims = [0 jend]+1;
myYtick = [0:dj:jend]+1;
%...............
zdepths = [0:dz:zzend];
for zj=1:length(zdepths)
    myYtickLabel{zj} = num2str(zdepths(zj));
end
%...............
set(ax1,'Ylim',myYlims,'Ytick',myYtick,'YTickLabel',myYtickLabel,'FontSize',SizeFontYticks)
set(ax1,'Xlim',myXlims,'XTick',myXtick,'XTickLabel',myXtickLabel,'Fontsize',SizeFontXticks)
%...............
% $$$ set(hcbar001,'Xticklabel',[]) %for 'horiz' colorbar.
% $$$ set(hcbar001,'Yticklabel',[]) %for 'vertic' colorbar.
%...............
set(hcbar001,'Visible','off'); %when using ax2 (otherwise, uncomment-line!!!!!!)
%...............
if strcmp(showlabelY,'not')
    set(ax1,'Yticklabel',[])
end
%...............
if strcmp(showlabelX,'not')
    set(ax1,'Xticklabel',[])
end
%...............
ax1TickLength=get(ax1,'ticklength');
ax1TickLengthNew=(1/8)*ax1TickLength;
set(ax1,'ticklength',ax1TickLengthNew);
%...............
%***************************************
% $$$ ax2=nan;
% $$$ hcbar002=nan;
% $$$ return

%%%%%%%%%%
%2nd AXES:
%%%%%%%%%%
ax2 = axes('Position', ax1position,'Visible', 'off');
%..............
% $$$ if minA~=maxA
% $$$     himg2=imagesc(A,[minA maxA]);
% $$$ else
% $$$     himg2=imagesc(A);
% $$$ end
%..............
if strcmp(keyPcolor,'not')
    %..............
    himg2=imagesc(A);
    %..............
    %%himg2=imagesc(Atmp);
    %..............
elseif strcmp(keyPcolor,'yes')
    himg2=pcolor(flipud(A));
    shading interp
end
set(gca,'clim',[Amin Amax]);
%..............
colormap(MAP) 
if strcmp(keyOctave,'yes')
    hcbar002 = colorbar(barpos);
else
    hcbar002 = mycolorbar(barpos);
end
%...............
myXtick=monthlim;
%...............
set(ax2,'Ylim',myYlims,'Ytick',myYtick,'YTickLabel',[])
set(ax2,'Xlim',myXlims,'XTick',myXtick,'XTickLabel',[])
% $$$ %%set(ax2,'ticklength',ax1TickLength);
% $$$ set(gca,'Fontsize',sizefont)
grid on
%...............
% $$$ ColorbarPosition=get(hcbar002,'Position');
% $$$ ColorbarPositionNew=ColorbarPosition;
% $$$ ColorbarPositionNew(1)=ColorbarPositionNew(1);
% $$$ ColorbarPositionNew(2)=ColorbarPositionNew(2);
% $$$ ColorbarPositionNew(3)=ColorbarPositionNew(3)/2;
% $$$ ColorbarPositionNew(4)=ColorbarPositionNew(4)/2;
% $$$ set(hcbar002,'Position',ColorbarPositionNew)
%..............
% $$$ if strcmp(keyOctave,'yes')
% $$$     if minA~=maxA
% $$$ 	set(colorbar(barpos),'Fontsize',sizefont,'Ylim',[minA maxA])
% $$$ 	set(colorbar(barpos),'Fontsize',sizefont,'Xlim',[minA maxA])
% $$$     else
% $$$ 	set(colorbar(barpos),'Fontsize',sizefont)
% $$$     end
% $$$ else
% $$$     if minA~=maxA
% $$$ 	set(mycolorbar(barpos),'Fontsize',sizefont,'Ylim',[minA maxA])
% $$$ 	set(mycolorbar(barpos),'Fontsize',sizefont,'Xlim',[minA maxA])
% $$$     else
% $$$ 	set(mycolorbar(barpos),'Fontsize',sizefont)
% $$$     end
% $$$ end
%..............
nanA=find(isnan(A(:))==1);
%..............
if maxA~=0 & maxA~=-999 & length(nanA)~=mm*nn %if there is data in mX.
    hclim=[minA maxA];
    if strcmp(ColorbarPosition,'Horizontal')
	set(hcbar002,'Xlim',hclim,'Fontsize',SizeFontColorbar)
    elseif strcmp(ColorbarPosition,'Vertical')
	set(hcbar002,'Ylim',hclim,'Fontsize',SizeFontColorbar)
    end
else %if there is no data in mX (ie. all -999).
    if strcmp(ColorbarPosition,'Horizontal')
	set(hcbar002,'Xlim',[0 eps],'Xtick',[],'Xticklabel',[],'Fontsize',SizeFontColorbar)
    elseif strcmp(ColorbarPosition,'Vertical')
	set(hcbar002,'Ylim',[0 eps],'Ytick',[],'Yticklabel',[],'Fontsize',SizeFontColorbar)
    end
end
%%set(hcbar002,'Fontsize',sizefont)
%..............
varargout{1}=himg2;
%..............

if strcmp(keyColorbar,'not') 
% $$$     delete(hcbar002) %(remove the colorbar)
    set(hcbar001,'Visible','off'); %(make invisible the colorbar)
    set(hcbar002,'Visible','off'); %(make invisible the colorbar)
end

%***************************
return

%..............................................
%ANTES LO HACIA ASI, PERO CREO QUE ES MAS SIMPLE DE ENTENDER COMO ESTA ARRIBA:
% $$$ himg1=imagesc(A);
% $$$ hc = mycolorbar(barpos);
% $$$ set(gca,'Xlim',[1,ndays],'XTick',monthlim2,'XTickLabel',xticklabels)
% $$$ set(gca,'Ylim',[0,zzend],'YTick',[0:dz:zzend],'Yticklabel',[0:dz:zzend])
% $$$ set(gca,'Fontsize',[5])
% $$$ set(hc,'Fontsize',[5])
% $$$ % $$$ axis square
% $$$ original = gca;
% $$$ position = get(original, 'Position');
% $$$ temporary = axes('Position', position,'Visible', 'off');
% $$$ himg2=imagesc(A);
% $$$ set(temporary,'Xlim',[1,ndays],'XTick',monthlim,'XTicklabel','');
% $$$ set(temporary,'Ylim',[0,zzend],'YTick',[0:dz:zzend],'Yticklabel','');
% $$$ % $$$ axis square
% $$$ grid on
%..............................................
