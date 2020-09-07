function [Mchg]=myChgResolutionIncrease(M,wi,wj)

%*********************************************************************
%MYCHGRESOL.m: Este programa aumenta la resolucion de X(m,n) a Xout(mm,nn), 
%donde "m<mm" y "n<nn" (o sea, aumento la resolucion: paso de pocos pixeles
%grandes a muchis pixeles pequenos). Se basa en repetir sobre grupos de
%pixeles grandes para obtener muchos pequenos iguales a cada grande, de modo que la
%size de X y Xout debe ser un multiplo entero.
%
%"wi" es el numero de pixels que entran en lat.
%"wj" es el numero de pixels que entran en log. 
%(de modo que "wi*wj" es el numero nuevo de pixels pequenos que repiten el grande)
%
% [Mchg]=myIncreaseResolution(M,wi,wj)
%
%*********************************************************************
[m,n,p]=size(M); %Low resolution initial.
I=[(1/wi):(1/wi):m]; %High resolution ouptut.
J=[(1/wj):(1/wj):n];
mm=length(I);
nn=length(J);

mnp=[m,n,p]
mmnnpp=[mm,nn,p]

if mod(mm,m)~=0 | mod(nn,n)~=0
    [m,mm]
    [n,nn]
    display('Ojo!, cambia el valor de "w"')
    error('tanto m y mm como n y nn deben ser multiplos enteros!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. PASO LOS DATOS DE (m,n) A UNA GRILLA DE (mm,nn):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:p %meses
    k
    ii=1;
    for i=1:m
	jj=1;
	for j=1:n
	    ij=[i,j];
	    %.................
	    X=M(i,j,k);
	    Xm=ones(wi,wj)*X;
	    %.................
	    wlat = ii:ii+(wi-1);
	    wlon = jj:jj+(wj-1);
	    %.................
	    Mchg(wlat,wlon,k)=Xm;
	    %.................
	    jj=jj+wj;
	end
	ii=ii+wi;
    end
end
