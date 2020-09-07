function [Xchg]=myChgResolutionDecrease(X,wi,wj)

%*********************************************************************
%MYCHGRESOL.m: Este programa reduce la resolucion de X(m,n) a Xout(mm,nn), 
%donde "m > mm" y "n > nn" (o sea, disminuyo la resolucion: paso de muchos pixeles
%pequenos a pocos pixeles grandes). Se basa en promediar sobre grupos de
%pixeles pequenos para obtener uno grande para cada grupo, de modo que la
%size de X y Xout debe ser un multiplo entero.
%
%"wi" es el numero de pixels que entran en lat.
%"wj" es el numero de pixels que entran en log. 
%(de modo que "wi*wj" es el numero de pixels que entran en la media para dar un pixel gordo). 
%
% [Xout]=mychgresol(X,wi,wj)
%
%*********************************************************************
[m,n,p]=size(X);
I=[1:wi:m];
J=[1:wj:n];
mm=length(I);
nn=length(J);

if mod(m,mm)~=0 | mod(n,nn)~=0
    [m,mm]
    [n,nn]
    display('Ojo!, cambia el valor de "w"')
    error('tanto m y mm como n y nn deben ser multiplos enteros!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. PASO LOS DATOS DE (m,n) A UNA GRILLA DE (mm,nn):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:p %time
    ktime = k 
    ii=0;
    for i=1:wi:m
	ii=ii+1;
	jj=0;
	%.................
% $$$ 	if mod(ii,90)==0
% $$$ 	    ilat = ii
% $$$ 	end
	%.................
	for j=1:wj:n
	    %.................
	    jj=jj+1;
	    [ii,jj];
	    %.................
	    wlat = i:i+(wi-1);
	    wlon = j:j+(wj-1);
	    %.................
	    Xw=X(wlat,wlon,k);
	    %.................
	    xave=nanmean(Xw(:));
	    xmed=nanmedian(Xw(:));
	    %.................
	    Xmean(ii,jj,k)=xave;
	    Xmedn(ii,jj,k)=xmed;
	    %.................
	end
    end
end

%%%%%%%%
%OUTPUT:
%%%%%%%%
Xchg = Xmean;
% $$$ Xchg = Xmedn;
