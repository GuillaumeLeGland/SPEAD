function [r2,a,b,nptos,pval]=myregress(y,x) %si quiero "sombra" con std(y).

%**********************************
%Syntaxis: [r2,corte,pente,npuntos,pval]=myregress(y,x)
%
% type='shadow': Shadow-std.
% type='errorbar': Errorbar-std.
%**********************************
y=y(:);
x=x(:);
[ychap,a,b,bs,aint,bint,res,resint,r2,F,pval,nptos]=sregress(y,x);

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
