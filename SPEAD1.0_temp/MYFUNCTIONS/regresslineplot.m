function [r2,corte,pente,npuntos]=regresslineplot(y,x)
%**********************************
%Syntaxis: [r2,corte,pente,npuntos]=regresslineplot(y,x)
%**********************************
y=y(:);
x=x(:);
[ychap,a,b,bs,aint,bint,res,resint,r2,F,pval,nptos]=sregress3(y,x);

%QUITO POSIBLES VALORES NEGATIVOS DE LA RECTA DE REGRESION:
J=find(ychap<0);
ychap(J)=nan;
x(J)=nan;
%INTERPOLO LA YCHAP (QUITO PRIMERO LOS POSIBLES NAN):
I=find(isnan(x)==0);
xx=x(I);
xxi=[min(xx):0.01:max(xx)];
J=find(isnan(x)==0 & isnan(ychap)==0);
yychap=ychap(J);
xx=x(J);
ychapI=interp1(xx,yychap,xxi);
%PLOTEO:
hp=plot(x,y,'r.',xxi,ychapI,'b-');
set(gca,'FontSize',[6]);
npuntos=num2str(nptos);
factor=100; %USAR ESTE NORMALMENTE (2 DECIMALES).
% $$$ factor=1000;
corte=num2str(round(factor*a)/factor);
if b>=0
    pente=['+ ',num2str(round(factor*b)/factor)];
else
    pente=['- ',num2str(abs(round(factor*b)/factor))];
end
r2=num2str(round(factor*r2)/factor);    
