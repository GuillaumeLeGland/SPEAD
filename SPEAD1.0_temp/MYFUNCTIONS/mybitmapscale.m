function [Xstar,zdata,zbitmap,zbitmapinf,zbitmapsup,zlims,Itick]=mybitmapscale(X,Xlims,keyLogBase,varargin)
%***********************************************
%Use: [Xstar,zbitmap]=mybitmapscale(X,xmin,xmax)
%***********************************************

%.................................
% $$$ keyLogBase='Log2';
% $$$ %.................................
% $$$ % $$$ keyLogBase='Log10';
% $$$ %.................................
% $$$ % $$$ keyLogBase,
% $$$ % $$$ pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QUITO POSIBLES CEROS Y PONGO NaN:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J=find(X==0);
X(J)=nan;
xminorg=min(X(:));
xmaxorg=max(X(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO MAPA DE BITS (0-256):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nvarargin=length(varargin);
if nvarargin==0 %dejo los valores min-max originales de los datos.
    xmin=xminorg;
    xmax=xmaxorg;
elseif nvarargin==1 %impongo yo el valor minimo (y dejo el max. original).
    xmax=xmaxorg;
    xmin=varargin{1};
    %..................................
% $$$     if xminorg>xmin %si el limite min. real es superior al impuesto, cojo el real.
% $$$ 	xmin=xminorg;
% $$$     end
    %..................................
    I=find(X<xmin);
    X(I)=xmin;
elseif nvarargin==2 %impongo yo el valor min. y max.
    xmin=varargin{1};
    xmax=varargin{2};
    %..................................
% $$$     if xminorg>xmin %si el limite min. real es superior al impuesto, cojo el real.
% $$$ 	xmin=xminorg;
% $$$     end
% $$$     if xmaxorg<xmax %si el limite max. real es inf. al impuesto, cojo el real.
% $$$ 	xmax=xmaxorg;
% $$$     end
    %..................................
    I=find(X<xmin);
    J=find(X>xmax);
    X(I)=xmin;
    X(J)=xmax;
end
xminmax=[xmin,xmax];
%.................................
if strcmp(keyLogBase,'Log10')
    a=log10(xmin);
    b=(log10(xmax)-a)/256;
    Xstar=(log10(X)-a)/b;
elseif strcmp(keyLogBase,'Log2')
    a=log2(xmin);
    b=(log2(xmax)-a)/256;
    Xstar=(log2(X)-a)/b;
end
%.................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINO LOS TICKS DE LA COLORBAR PARA LOS BIT-MAPS (0-256):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%........................
xxmin=ceil(xmin*100)/100;
zdata=[xxmin:xmax];
%........................
if strcmp(keyLogBase,'Log10')
    %........................
    x1=0.1;
    x2=0.2;
    x3=0.3;
    x4=0.4;
    x5=0.5;
    %........................
    x=[x1,x3]; %Simetric spacing between ticks [0.1, 0.3, 1.0, 3.0, 10, 30, etc]
% $$$     x=[x2,x3];
% $$$     x=[x1,x2,x3];
% $$$     x=[x1,x2,x3,x4]; %For PreySwitch paper (Units PO4)
% $$$     x=[x1,x2,x3,x4,x5];
    %........................
    antilog10x=10.^(x);
    %........................
    z0000 = log10(antilog10x)*0.0001;
    z000  = log10(antilog10x)*0.001;
    z00   = log10(antilog10x)*0.01;
    z0    = log10(antilog10x)*0.1;
    z1    = log10(antilog10x)*1;
    z2    = log10(antilog10x)*10;
    z3    = log10(antilog10x)*100;
    z4    = log10(antilog10x)*1000;
    z5    = log10(antilog10x)*10000;
    %........................
    zdata=[z0000,z000,z00,z0,z1,z2,z3,z4,z5];
    %........................
    %Como tengo problemas de precision (no consigo por ej. hacer
    %i=find(z0==0.01)), debo hacer lo siguiente para tener decimales que se
    %puedan buscar y encontrar:
    zz0000 = (round(z0000*100000))/100000;
    zz000  = (round(z000*10000))/10000;
    zz00   = (round(z00*1000))/1000;
    zz0    = (round(z0*1000))/1000;
    zz1    = (round(z1*1000))/1000;
    zz2    = (round(z2*1000))/1000;
    zz3    = (round(z3*1000))/1000;
    zz4    = (round(z4*1000))/1000;
    zz5    = (round(z5*1000))/1000;
    zdata=[zz0000,zz000,zz00,zz0,zz1,zz2,zz3,zz4,zz5];
elseif strcmp(keyLogBase,'Log2')
    %........................
    x1=0.1;
    x2=0.2;
    x3=0.3;
    x4=0.4;
    x5=0.5;
    %........................
    x=[x1];
% $$$     x=[x1,x2,x3];
% $$$     x=[x1,x2,x3,x4,x5];
    %........................
    antilog2x = 2.^(x);
    %........................
    z6B    = log2(antilog2x)*10/64;
    z5B    = log2(antilog2x)*10/32;
    z4B    = log2(antilog2x)*10/16;
    z3B    = log2(antilog2x)*10/8;
    z2B    = log2(antilog2x)*10/4;
    z1B    = log2(antilog2x)*10/2;
    z0     = log2(antilog2x)*10;
    z1A    = log2(antilog2x)*10*2;
    z2A    = log2(antilog2x)*10*4;
    z3A    = log2(antilog2x)*10*8;
    z4A    = log2(antilog2x)*10*16;
    z5A    = log2(antilog2x)*10*32;
    z6A    = log2(antilog2x)*10*64;
    %........................
    zdata=[z6B,z5B,z4B,z3B,z2B,z1B,z0,z1A,z2A,z3A,z4A,z5A,z6A];
    %........................
end
zdataPto0 = zdata

%DEFINO DE DONDE-A-DONDE DEBE IR 'ZDATA':
%======================================
% $$$ %......................
% $$$ % $$$ %J=find(zdata>=xxmin & zdata<=round(xmax));
% $$$ % $$$ J=find(zdata>=xxmin & zdata<=ceil(xmax));
% $$$ J=find(zdata>=xmin & zdata<=ceil(xmax));
% $$$ %......................
% $$$ J1=J(1)-1;
% $$$ %J1=[J(1)-2,J(1)-1];
% $$$ J2=J(end)+1;
% $$$ %......................
% $$$ if J1~=0
% $$$     JJ=[J1,J,J2]; %cojo tb posic. justo anterior y posterior.
% $$$ elseif J1==0
% $$$     JJ=[J,J2];
% $$$ end
% $$$ %......................
% $$$ %zdata=zdata(J);
% $$$ zdata=zdata(JJ)
% $$$ %......................
%======================================
zdataPto1 = zdata
%......................
%Obligo a empezar 'zdata' por el xmin:
% $$$ I=find(zdata<=xmin); %posiciones en 'zdata' con valores inferiores (o igual) a xmin.
% $$$ if isempty(I)==0 %si no esta vacia.
% $$$     i=I(end); %posicion justo anterior a xmin (o pos. de xmin).
% $$$     zdata=zdata(i+1:end); %posicion justo en xmin.
% $$$ end
%......................
%Obligo a empezar 'zdata' por el xmin o un pelin por debajo:
I = find(zdata <= xmin); %posiciones en 'zdata' con valores inferiores (o igual) a xmin.
if isempty(I) == 0 %si no esta vacia.
    i = I(end); %posicion justo anterior a xmin (o pos. de xmin).
    zdata = zdata(i:end);
end
%......................
% $$$ %Obligo a empezar 'zdata' casi por el xmin:
% $$$ I=find(zdata<=xmin); %posiciones en 'zdata' con valores inferiores (o igual) a xmin.
% $$$ if isempty(I)==0 %si no esta vacia.
% $$$     i=I(end-1); %posicion dos veces anterior a xmin.
% $$$     zdata=zdata(i:end);
% $$$ end
%......................
%======================================
zdataPto2 = zdata
%......................
% $$$ %Obligo a acabar 'zdata' por el xmax:
% $$$ J=find(zdata>=xmax); %posiciones en 'zdata' con valores sup (o igual) a xmax.
% $$$ if isempty(J)==0
% $$$     j=J(1); %posicion justo posterior a xmax (o pos. de xmax).
% $$$     zdata=zdata(1:j-1); %posicion just en xmax.
% $$$ end
%......................
%Obligo a acabar 'zdata' por el xmax o un pelin por encima:
J=find(zdata>=xmax); %posiciones en 'zdata' con valores sup (o igual) a xmax.
if isempty(J)==0
    j=J(1); %posicion justo posterior a xmax (o pos. de xmax).
    zdata=zdata(1:j);
end
%......................
% $$$ %Obligo a acabar 'zdata' por un pelin encima del xmax:
% $$$ J=find(zdata>=xmax); %posiciones en 'zdata' con valores sup (o igual) a xmax.
% $$$ if isempty(J)==0
% $$$     j=J(2); %posicion dos veces posterior a xmax.
% $$$     zdata=zdata(1:j);
% $$$ end
%......................
%======================================
zdataPto3 = zdata
%......................
xminmax=[xmin,xmax] %show this.
%......................
pause(1)

%%%%%%%%%%%%%%%%%%%%
%DEFINE TICK LABELS:
%%%%%%%%%%%%%%%%%%%%
%======================
%......................
ztick=zdata;
%......................
nsize = length(ztick);
%......................
%======================
% $$$ %......................
% $$$ Inan = [1:2:nsize];
% $$$ %......................
% $$$ ztick(Inan)=nan;
% $$$ %......................
% $$$ Itick = find(isnan(ztick)==0);
% $$$ %......................
%======================
%......................
Itick = 1; %First tick ON!
% $$$ Itick = []; %First tick OFF!
%......................
dtick = 1; %One ticklabel every one tick.
% $$$ dtick = 2; %One ticklabel every two ticks.
%......................
for j=1:nsize
    if j==1
% $$$ 	index = 0; %For ASLOasm2012.
	index = 1; %For PreySwitch.
    end
    if dtick == 1
	index = index + 1;
    elseif dtick == 2
	if mod(j,dtick)==1; %Odd numbers
	    index = index + 2;
	else mod(j,dtick)==0; %Even numbers
	    index = index + 3;
	end
    end
    if index > nsize
	break %stop the for loop.
    end
    Itick = [Itick,index];
end
%======================
%......................
Itick
%......................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO LOS VALORES CON UNIDADES QUE CORRESPONDEN A CADA TICK DEL BIT-MAP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%================================
%................................
if strcmp(keyLogBase,'Log10')
    zbitmap    = (log10(zdata)-a)/b;
    zbitmapinf = (log10(xmin )-a)/b;
    zbitmapsup = (log10(xmax )-a)/b;
    zlims      = (log10(Xlims)-a)/b;
elseif strcmp(keyLogBase,'Log2')
    zbitmap    = (log2(zdata)-a)/b;
    zbitmapinf = (log2(xmin )-a)/b;
    zbitmapsup = (log2(xmax )-a)/b;
    zlims      = (log2(Xlims)-a)/b;
end
%................................
%================================
% $$$ %................................
% $$$ if strcmp(keyLogBase,'Log10')
% $$$     zbitmap=(log10(ztick)-a)/b;
% $$$     zbitmapinf=(log10(xmin)-a)/b;
% $$$     zbitmapsup=(log10(xmax)-a)/b;
% $$$ elseif strcmp(keyLogBase,'Log2')
% $$$     zbitmap=(log2(ztick)-a)/b;
% $$$     zbitmapinf=(log2(xmin)-a)/b;
% $$$     zbitmapsup=(log2(xmax)-a)/b;
% $$$ end
% $$$ %................................
% $$$ Inonan = find(isnan(ztick)==0);
% $$$ %................................
% $$$ zdata=zdata(Inonan)
% $$$ zbitmap=zbitmap(Inonan)
% $$$ %................................
%================================
ZandZbitmap=[zdata(:),zbitmap(:)] %show this.

%**********************************
% $$$ display('END mybitmapscale.m')

%....
% $$$ figure(100)
% $$$ MAP=jet;
% $$$ MAP(1,:)=[0 0 0];
% $$$ 
% $$$ subplot(1,2,1)
% $$$ imagesc(X)
% $$$ colormap(MAP)
% $$$ hc=colorbar('horiz');
% $$$ axis square
% $$$ 
% $$$ subplot(1,2,2)
% $$$ imagesc(Xstar)
% $$$ colormap(MAP)
% $$$ hc=colorbar('horiz');
% $$$ %set(hc,'Xlim',[zbitmap(1) zbitmap(end)],'Xtick',zbitmap(:),'Xticklabel',zdata)
% $$$ set(hc,'Xlim',[zbitmapinf zbitmapsup],'Xtick',zbitmap(:),'Xticklabel',zdata) %USAR ESTE!!!!
% $$$ axis square
