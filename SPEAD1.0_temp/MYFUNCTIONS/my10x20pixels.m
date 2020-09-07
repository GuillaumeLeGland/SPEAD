function [Avar,Cavar]=my10x20pixels(deltalat,deltalon,VAR,varargin)

%******************************
%Program MY10X20PIXELS.m:
%Este programa obtiene una matriz 3D (18,18,12) a partir de una matriz
%3D(180,360,12) por medio de promediar los datos en pixels grandes de (10 x 20).
%La variable input (VAR) no necesita obligatoriamente ser una matriz 3D;
%tambien puede ser una columna de datos, pero entonces se debe incluir
%como inputs 3 columnas mas con el MES, LATITUD, LONGITUD de cada dato.
%
%Use: [Avar,Cavar]=dmsalgoritmo_10x20pixels(deltalan,deltalon,VAR,MESES,LATITUD,LONGITUD)
%******************************
nvarargin=length(varargin);
if nvarargin>0
    MESES=varargin{2};
    LATITUD=varargin{3};
    LONGITUD=varargin{4};
end

[m,n,p]=size(VAR);

%...........
% $$$ deltalat=10;
% $$$ deltalon=20;
%...........
% $$$ deltalat=5;
% $$$ deltalon=10;
%...........
LATRG=[+90:-deltalat:-90];
LONRG=[-180:+deltalon:+180];
mm=length(LATRG)-1; %numero de intervalos de latitud (son 18).
nn=length(LONRG)-1; %numero de intervalos de longitud (son 18).
for k=1:12
    k
    if p==1 %matriz 2D.
	K=find(MESES==k);
	LATITUDk=LATITUD(K);
	LONGITUDk=LONGITUD(K);
	VARk=VAR(K);
	for i=1:mm
	    latsup=LATRG(i);
	    latinf=LATRG(i+1);
	    for j=1:nn
		lonsup=LONRG(j+1);
		loninf=LONRG(j);
		J=find(LATITUDk>latinf & LATITUDk<=latsup & LONGITUDk>loninf & LONGITUDk<=lonsup);
		VARkJ=VARk(J);
		var=nanmean(VARkJ); %VAR promedio en el (10x20)size pixel para mesk.
		npvar=length(find(VARkJ>=0));
		if npvar>0
		    c=npvar; %Numero de puntos usados en el promedio.
		else
		    c=nan;
		end
		Avar(i,j,k)=var;
		Cavar(i,j,k)=c;
	    end
	end
    elseif p==12 %matriz 3D.
	VARk=VAR(:,:,k);
	for i=1:mm
	    imin=(i*deltalat)-(deltalat-1);
	    imax=(i*deltalat);
	    for j=1:nn
		jmin=(j*deltalon)-(deltalon-1);
		jmax=(j*deltalon);
		VARkpixel=VARk(imin:imax,jmin:jmax);
 		var=nanmean(VARkpixel(:)); %USAR ESTE!!!!
		npvar=length(find(VARkpixel(:)>0));
		if npvar>0
		    c=npvar; %Numero de puntos usados en el promedio.
		else
		    c=nan;
		end
		Avar(i,j,k)=var;
		Cavar(i,j,k)=c;
	    end
	end
    end
end
