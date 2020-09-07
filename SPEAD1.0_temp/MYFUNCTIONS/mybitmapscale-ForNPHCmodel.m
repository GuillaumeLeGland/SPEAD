function [Xstar,z,zbitmap,zbitmapinf,zbitmapsup]=mybitmapscale(X,varargin)
%***********************************************
%Use: [Xstar,zbitmap]=mybitmapscale(X,xmin,xmax)
%***********************************************
    
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
a=log10(xmin);
b=(log10(xmax)-a)/256;
Xstar=(log10(X)-a)/b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINO LOS TICKS DE LA COLORBAR PARA LOS BIT-MAPS (0-256):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xxmin=ceil(xmin*100)/100;
z=[xxmin:xmax];
x1=0.1;
x2=0.2;
x3=0.3;
x4=0.4;
x5=0.5;
x=[x1,x2,x3,x4,x5];
antilog10x=10.^(x);
z00=log10(antilog10x)*0.01;
z0=log10(antilog10x)*0.1;
z1=log10(antilog10x)*1;
z2=log10(antilog10x)*10;
z3=log10(antilog10x)*100;
z4=log10(antilog10x)*1000;
z5=log10(antilog10x)*10000;
z=[z00,z0,z1,z2,z3,z4,z5];
%Como tengo problemas de precision (no consigo por ej. hacer
%i=find(z0==0.01)), debo hacer lo siguiente para tener decimales que se
%puedan buscar y encontrar:
zz00=(round(z00*1000))/1000;
zz0=(round(z0*1000))/1000;
zz1=(round(z1*1000))/1000;
zz2=(round(z2*1000))/1000;
zz3=(round(z3*1000))/1000;
zz4=(round(z4*1000))/1000;
zz5=(round(z5*1000))/1000;
z=[zz00,zz0,zz1,zz2,zz3,zz4,zz5];

%DEFINO DE DONDE-A-DONDE DEBE IR 'Z':
%......................
% $$$ % $$$ %J=find(z>=xxmin & z<=round(xmax));
% $$$ % $$$ J=find(z>=xxmin & z<=ceil(xmax));
% $$$ J=find(z>=xmin & z<=ceil(xmax));
% $$$ J1=J(1)-1;
% $$$ %J1=[J(1)-2,J(1)-1];
% $$$ J2=J(end)+1;
% $$$ if J1~=0
% $$$     JJ=[J1,J,J2]; %cojo tb posic. justo anterior y posterior.
% $$$ elseif J1==0
% $$$     JJ=[J,J2];
% $$$ end
% $$$ %z=z(J);
% $$$ z=z(JJ)
%......................
%======================================
%......................
%Obligo a empezar 'z' por el xmin:
I=find(z<=xmin); %posiciones en 'z' con valores inferiores (o igual) a xmin.
if isempty(I)==0 %si no esta vacia.
    i=I(end); %posicion justo anterior a xmin (o pos. de xmin).
    z=z(i:end);
end
%......................
%Obligo a empezar 'z' casi por el xmin:
% $$$ I=find(z<=xmin); %posiciones en 'z' con valores inferiores (o igual) a xmin.
% $$$ if isempty(I)==0 %si no esta vacia.
% $$$     i=I(end-1); %posicion dos veces anterior a xmin.
% $$$     z=z(i:end);
% $$$ end
%......................
%======================================
%......................
%Obligo a acabar 'z' por el xmax:
J=find(z>=xmax); %posiciones en 'z' con valores sup (o igual) a xmax.
if isempty(J)==0
    j=J(1); %posicion justo posterior a xmax (o pos. de xmax).
    z=z(1:j);
end
%......................
%Obligo a acabar 'z' por un pelin encima del xmax:
% $$$ J=find(z>=xmax); %posiciones en 'z' con valores sup (o igual) a xmax.
% $$$ if isempty(J)==0
% $$$     j=J(2); %posicion dos veces posterior a xmax.
% $$$     z=z(1:j);
% $$$ end
%......................
%======================================
xminmax=[xmin,xmax]; %show this.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO LOS VALORES CON UNIDADES QUE CORRESPONDEN A CADA TICK DEL BIT-MAP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zbitmap=(log10(z)-a)/b;
zbitmapinf=(log10(xmin)-a)/b;
zbitmapsup=(log10(xmax)-a)/b;
ZandZbitmap=[z(:),zbitmap(:)]; %show this.

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
% $$$ %set(hc,'Xlim',[zbitmap(1) zbitmap(end)],'Xtick',zbitmap(:),'Xticklabel',z)
% $$$ set(hc,'Xlim',[zbitmapinf zbitmapsup],'Xtick',zbitmap(:),'Xticklabel',z) %USAR ESTE!!!!
% $$$ axis square
