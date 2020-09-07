function [ts]=fishertest(ro,r,n)

%Test de FISHER para evaluar si 2 valores de correlacion son diferentes:
%[ts]=fishertest(ro,r,n)
%Ho: ro = ro
%H1: r ~= ro.
%n = number of data.
%(see SOAKEL, pag. 585.)
    
%Z-TRANSFORMATION OF "r":
%--------------------------------------------
%z = 1/2 * Ln * ((1+r)/(1-r))
%if r=0.837 => z=1.2111
%if r=0.500 => z=0.5493
%--------------------------------------------
fi=(1/2)*log((1+ro)/(1-ro));
z=(1/2)*log((1+r)/(1-r));
[fi,z]

%t-TEST:
%--------------------------------------------
%ts = (z-fi)/(1/sqrt(n-3)) = (z-fi)*sqrt(n-3)
%--------------------------------------------
ts = (z-fi)*sqrt(n-3);


















