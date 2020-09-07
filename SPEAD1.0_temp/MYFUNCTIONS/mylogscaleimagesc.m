function [hbar]=mylogscaleimagesc(X,fignum)

%*******************************************
%Program MYLOGSCALEIMAGESC.m:
%This programs takes a 2D(m,n) array of data and makes an "imagesc"
%figure of the data in log10 scale, with its corresponding "colorbar"
%using ticklabels with original units.
%
%Use: mylogscaleimagesc(X,'fignum')
%
% X: 2D(m,n) array of original data.
% fignum: number of the figure (use '1' for example).
%*******************************************

%colorbartype='horiz';
colorbartype='vertic';

%%%%%%%%%%%%%%%%%%%%%
%CHANGE ZEROS BY NAN:
%%%%%%%%%%%%%%%%%%%%%
%....................
pcnt = 0.05;
%....................
xmean = nanmean(X(:));
%....................
% $$$ J=find(X==0);
%....................
J=find(X <= pcnt*xmean);
%....................
X(J)=nan;
xmin=min(X(:));
xmax=max(X(:));
minmaxX=[xmin,xmax];

%%%%%%%%%%%%%%%%%%%%%%
%MAKE BIT-MAP (0-256):
%%%%%%%%%%%%%%%%%%%%%%
a=log10(xmin);
b=(log10(xmax)-a)/256;
Xstar=(log10(X)-a)/b;

%%%%%%%%%%%%%%%%%%
%VARIABLE = Xstar:
%%%%%%%%%%%%%%%%%%
Array2D=Xstar;
amax=max(Array2D(:));
amin=min(Array2D(:));
aminmax=[amin,amax];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SET DE 'TICKS' FOR THE COLORBAR OF THE BIT-MAP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%..........................
xxmin=ceil(xmin*100)/100;
z=[xxmin:xmax];
%..........................
x1=0.1;
x2=0.2;
x3=0.3;
x4=0.4;
x5=0.5;
%..........................
x=[x1,x2,x3,x4,x5];
antilog10x=10.^(x);
%..........................
z00=log10(antilog10x)*0.01;
z0=log10(antilog10x)*0.1;
z1=log10(antilog10x)*1;
z2=log10(antilog10x)*10;
z3=log10(antilog10x)*100;
z4=log10(antilog10x)*1000;
z5=log10(antilog10x)*10000;
z6=log10(antilog10x)*100000;
%..........................
%To avoid precision problems:
zz00=(round(z00*1000))/1000;
zz0=(round(z0*1000))/1000;
zz1=(round(z1*1000))/1000;
zz2=(round(z2*1000))/1000;
zz3=(round(z3*1000))/1000;
zz4=(round(z4*1000))/1000;
zz5=(round(z5*1000))/1000;
zz6=(round(z6*1000))/1000;
%..........................
z=[zz00,zz0,zz1,zz2,zz3,zz4,zz5,zz6];
%..........................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SET MAXIMUM AND MINIMUM VALUE FOR 'Z':
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Force 'z' to start at xmin:
I=find(z<=xmin); %posic. in 'z' with values less (or equal) than xmin.
if isempty(I)==0 %if not empty.
    i=I(end); %posic. just before xmin (or pos. of xmin).
    z=z(i:end);
end

%Force 'z' to finish at xmax:
J=find(z>=xmax); %posic. in 'z' with values sup (or equal) than xmax.
if isempty(J)==0 %if not empty.
    j=J(1); %posic. just after xmax (or pos. of xmax).
    z=z(1:j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET THE REAL VALUES CORRESPONDING TO EACH TICK OF THE 'Z' BIT-MAP COLORBAR:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zbitmap=(log10(z)-a)/b;
zbitmapinf=(log10(xmin)-a)/b;
zbitmapsup=(log10(xmax)-a)/b;
zandzbitmap=[z;zbitmap];

%%%%%%%%%
%IMAGESC:
%%%%%%%%%
cmap=jet;
nancolor = [1 1 1];

%size of colormap:
nbits = size(cmap,1);

%color step:
dmap=(amax-amin)/nbits;

%add nan color to colormap:
MAP = [nancolor; cmap];
%%cmap(1,:)=[1 1 1];
%%cmap(1,:)=[0 0 0];
%%cmap(1,:)=[0.5 0.5 0.5];

figure(fignum)
himg=imagesc(Array2D);
hbar=colorbar(colorbartype);
colormap(MAP)

%changing color limits:
caxis([amin-dmap amax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE COLORBAR TICK AND TICKLABELS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(colorbartype,'horiz')
    set(hbar,'Xlim',[zbitmap(1) zbitmap(end)],'Xtick',zbitmap(:),'Xticklabel',z)
elseif strcmp(colorbartype,'vertic')
    set(hbar,'Ylim',[zbitmap(1) zbitmap(end)],'Ytick',zbitmap(:),'Yticklabel',z)
end

%change Y limit for colorbar to avoid showing NaN color
ylim(hbar,[amin amax])



