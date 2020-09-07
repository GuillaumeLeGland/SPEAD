function [Xpct]=mypercentile(X,pct) %NO USAR!!!! (usar "percentile.m")
%-----------------------------------------
%Use: [Xpct]=mypercentile(X,pct)
%pct: debe ser un numero entre 1 y 20.
%-----------------------------------------
%NO estoy seguro de que esto funcione correctamente, me lo baje de:
%<http://matlabdatamining.blogspot.com/search/label/percentiles>

x = X(:);
I = find(x >= 0);
x = x(I);

n = 20 %si quiero percentiles de 5% (100/20 = 5).
F = ceil(n*tiedrank(x)/length(x)); %divide X en 20 subgrupos de 5%.
tabulate(F)

% $$$ xi=x(F==pct); %
% $$$ Xpct=xi(1)

xi = x(F==(pct-1)); %
Xpct = xi(end)

