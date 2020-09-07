function [varargout]=myregresslineloglog(y,x,varargin) %si quiero "sombra" con std(y).
%function [r2,corte,pente,npuntos]=myregresslineplot(y,x,varargin) %si quiero "sombra" con std(y).
%**********************************
%Syntaxis: [r2,corte,pente,npuntos,pval]=myregresslineplot(y,x,stdy,type)
%
% type='shadow': Shadow-std.
% type='errorbar': Errorbar-std.
%**********************************
%..........
yorig=y(:);
xorig=x(:);
%..........
y=log10(yorig);
x=log10(xorig);
%..........

[ychap,a,b,bs,aint,bint,res,resint,r2,F,pval,nptos]=sregress(y,x);
%........................
if length(varargin)==1
    type='errorbar';
    stdy=varargin{1};
    stdy=stdy(:);
elseif length(varargin)==2
    type=varargin{2}
    stdy=varargin{1};
    stdy=stdy(:);
else
    type='NO';
end
%........................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULO EL STANDARD ERROR (SE) DE LA PENDIENTE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%<http://stattrek.com/AP-Statistics-4/Test-Slope.aspx?Tutorial=AP>
disty=y-ychap; %distancia de los puntos y a la recta de regress.
dispx=x-mean(x); %dispersion de los puntos x (distancia a su media).

sumdisty2=sum(disty.^2);
sumdispx2=sum(dispx.^2);

SEslope=sqrt(sumdisty2/(nptos-2))/sqrt(sumdispx2);

%%%%%%%%%%%%%%%%%%
%t-TEST STATISTIC:
%%%%%%%%%%%%%%%%%%
Ts=b/SEslope;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the p-values:
%%%%%%%%%%%%%%%%%%%%%%%%%
%....
%NOTA: Falla si hay NaNs...
% $$$ pvalslope=2*(1-tcdf(abs(Ts),nptos-2));
% $$$ 
% $$$ distpval=pval-pvalslope;
% $$$ if distpval>1d-10
% $$$     [SEslope,Ts,pvalslope,pval],pause
% $$$     display('error: mi calculo del pval de la slope NO sale igual que el del matlab!!')
% $$$     pause
% $$$ end
%....

%QUITO POSIBLES VALORES NEGATIVOS DE LA RECTA DE REGRESION:
% $$$ ychapold=ychap;
% $$$ J=find(ychap<0);
% $$$ ychap(J)=nan;
% $$$ x(J)=nan;
% $$$ [ychapold(:),ychap(:)],pause

%INTERPOLO LA YCHAP (QUITO PRIMERO LOS POSIBLES NAN):
I=find(isnan(x)==0);
xx=x(I);
xxmin=min(xx);
xxmax=max(xx);
%[xxmin,xxmax]
dx=(xxmax-xxmin)/100;
xxi=[xxmin:dx:xxmax];
J=find(isnan(x)==0 & isnan(ychap)==0);
yychap=ychap(J);
xx=x(J);
ychapI=interp1(xx,yychap,xxi);
%PLOTEO:
%....................
%SI QUIERO SOMBRA:
if strcmp(type,'shadow')
    Inonan=find(isnan(y)==0);
    xx=x(Inonan);
    yy=y(Inonan);
    sstdy=stdy(Inonan);
    [xx(:),yy(:),sstdy(:)]
    confplot(xx,yy,sstdy);
    hold on
elseif strcmp(type,'errorbar')
    %he=errorbar(x,y,stdy);
    he=v6errorbar(x,y,stdy);
    set(he(2),'LineStyle','none')
    hold on
end
%....................
% $$$ if strcmp(type,'shadow')
% $$$     confplot(x,y,stdy);
% $$$     hold on
% $$$ end
%....................
hp=loglog(xorig,yorig,'r.',10.^xxi,10.^(ychapI),'b-');
hold off
%....................
% $$$ if strcmp(type,'errorbar')
% $$$     hold on
% $$$     he=errorbar(x,y,stdy);
% $$$     hold off
% $$$     set(he(2),'LineStyle','none')
% $$$ end
%....................
% $$$ set(hp(1),'Markersize',[5]);
% $$$ set(hp(2),'LineWidth',[1]);
%....................
% $$$ set(hp(1),'Markersize',[15]);
% $$$ set(hp(2),'LineWidth',[1]);
%....................
%PARA SCIENCE PAPER:
set(hp(1),'Markersize',[20]);
set(hp(2),'Color',[0 0 0],'LineWidth',[1]);
%....................
%PARA DMS2006 MEETING:
% $$$ set(hp(1),'Markersize',[40]);
% $$$ set(hp(2),'Color',[0 0 0.5],'LineWidth',[2]);
%....................
%set(hp(1),'Marker','o',);
%set(hp(2),'Color',[0 0 0.5],'LineWidth',[2]);
%set(hp(2),'Color',[0 0.5 0],'LineWidth',[2]);
% $$$ set(hp(2),'Color',[0 0 0],'LineWidth',[2]);
%....................
% $$$ set(gca,'FontSize',[6]);
%....................
npuntos=num2str(nptos);
%factor=100; %USAR ESTE NORMALMENTE (2 DECIMALES).
factor=1000; %(3 DECIMALES).
corte=num2str(round(factor*a)/factor);
if b>=0
    pente=['+ ',num2str(round(factor*b)/factor)];
else
    pente=['- ',num2str(abs(round(factor*b)/factor))];
end
r2out=num2str(round(100*r2)/100);    
pvalout=num2str(round(1000*pval)/1000);    

%OUTPUT:
% $$$ varargout{1}=r2;
% $$$ varargout{2}=a;
% $$$ varargout{3}=b;
% $$$ varargout{4}=nptos;
% $$$ varargout{5}=pval;

%OUTPUTbis (EN STRING):
varargout{1}=r2out;
varargout{2}=corte;
varargout{3}=pente;
varargout{4}=npuntos;
varargout{5}=pvalout;
