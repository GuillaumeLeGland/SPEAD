function [VARI]=myPlankTOMregrid(VAR,Land)

%CONDICIONES DE FRONTERA:
for k=1:12
    VARk=VAR(:,:,k);
    %Anado 2 filas mas (son las condiciones de frontera
    %para la interpolacion): la repito a lat=0 y la
    %primera la repito a lat=181.
    W(:,:,k)=[VARk(:,end),VARk,VARk(:,1)]; 
end

%REGRID:
[m,n,p]=size(W);
%lat:
y=[1:152];
dy=152/180;
yi=[dy:dy:152];
%long:
x=[1:n];
dx=n/(2*n); %aumento al doble el numero de columnas (paso de 180+2 a 360+4).
xi=[dx:dx:n];
[X,Y]=meshgrid(x,y);
[XI,YI]=meshgrid(xi,yi);
for k=1:12
    k
    Wk=W(:,:,k);
    WIk = griddata(X,Y,Wk,XI,YI);
    WI(:,:,k)=WIk(:,1+2:end-2); %quito las columnas que use de frontera.
end

%CENTRO EN EL ATLANTICO:
for k=1:12
    WIk=WI(:,:,k);
    longlim=100;
    WIk=[WIk(:,longlim+1:end),WIk(:,1:longlim)];
    WIk(Land)=nan;
    VARI(:,:,k)=WIk;
end
