function [MOUT]=myfunchangeresol(M,w)

%*********************************************************************
%Este programa convierte la resolucion de X(m,n) a Xout(mm,nn), donde
%m>mm y n>nn (o sea, disminuyo la resolucion: paso de muchos pixeles
%pequenos a pocos pixeles grandes). Se basa en promediar sobre grupos de
%pixeles pequenos para obtener uno grande para cada grupo, de modo que la
%size de X y Xout debe ser un multiplo entero.
%"w" es el numero de pixels que entran en lat. y en log (de modo que
%"w*w" es el numero de pixels que entran en la media para dar un pixel gordo). 
%*********************************************************************
[m,n,p]=size(M); 
I=[1:w:m];
J=[1:w:n];
mm=length(I);
nn=length(J);

if mod(m,mm)~=0 | mod(n,nn)~=0
    display('Ojo!, cambia el valor de "w"')
    error('tanto m y mm como n y nn deben ser multiplos enteros!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. PASO LOS DATOS DE (m,n) A UNA GRILLA DE (mm,nn):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:p %meses
    ii=0;
    for i=1:w:m
	ii=ii+1;
	jj=0;
	for j=1:w:n
	    jj=jj+1;
	    X=M(i:i+(w-1),j:j+(w-1),k);
	    xm=nanmean(nanmean(X));
	    Xm(ii,jj,k)=xm;
	end
    end
end
MOUT=Xm;
