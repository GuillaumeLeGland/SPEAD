% $$$ function [xr,yr]=mypolyregress(x,y,dx,xlim,order);
function [yr]=mypolyregress(x,y,xr,order);
    
x=x(:);
y=y(:);

%**************************************************
% $$$ %Use: [xr,yr]=mypolyregress(x,y,dx,xlim,order);
%Use: [yr]=mypolyregress(x,y,xr,order);

%x: variable independiente.
%y: variable dependiente (es sobre la que quiero hacer la regres).
%dx: grid space.
% $$$ %xlim: grid max value.
%xr: profs. donde quiero obtener la curva de regresion.
%order: Polinomio de orden tal.
%**************************************************

%............................
% $$$ PolyRegressOrder=order
%............................
% $$$ xr=[1:dx:xlim]'; %profs. donde quiero obtener la curva de regresion.
%............................
xr=xr(:);
%............................

if order==1
    X=[ones(length(x),1),x];
    Xr=[ones(length(xr),1),xr];
elseif order==2
    X=[ones(length(x),1),x,x.^2];
    Xr=[ones(length(xr),1),xr,xr.^2];
elseif order==3
    X=[ones(length(x),1),x,x.^2,x.^3];
    Xr=[ones(length(xr),1),xr,xr.^2,xr.^3];
elseif order==4
    X=[ones(length(x),1),x,x.^2,x.^3,x.^4];
    Xr=[ones(length(xr),1),xr,xr.^2,xr.^3,xr.^4];
elseif order==5
    X=[ones(length(x),1),x,x.^2,x.^3,x.^4,x.^5];
    Xr=[ones(length(xr),1),xr,xr.^2,xr.^3,xr.^4,xr.^5];
elseif order==6
    X=[ones(length(x),1),x,x.^2,x.^3,x.^4,x.^5,x.^6];
    Xr=[ones(length(xr),1),xr,xr.^2,xr.^3,xr.^4,xr.^5,xr.^6];
end
b=(X'*X)\(X'*y); %coefficientes de regression.
yr=Xr*b;         %polynomio de regression.

%*******************
%-------------------------------
%Example: Hay un error que no entiendo:
% $$$ xlim=50;
% $$$ xmax=50;
% $$$ x=[[1:10],xmax];
% $$$ y=rand(1,11);
% $$$ [xr,yr]=mypolyregress(x(:),y(:),1,xlim,6);
% $$$ figure(1)
% $$$ plot(x,y,'r*',xr,yr,'-')
% $$$ axis([0 xmax, -inf +inf])
% $$$ 
% $$$ xmax=20;
% $$$ x=[[1:10],xmax];
% $$$ [xr,yr]=mypolyregress(x(:),y(:),1,xlim,6);
% $$$ figure(10)
% $$$ plot(x,y,'r*',xr,yr,'-')
% $$$ axis([0 xmax, -inf +inf])
%-------------------------------
