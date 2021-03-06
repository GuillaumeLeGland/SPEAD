function [xgm] = mygeomean1D(x)

%============================================================
%............................................................
if any(x < 0)
    error('The data must all be non-negative numbers.')
end
%............................................................
if all(x) %if all non-zero values.
    xs = x; 
else %if there are some zero values.
% $$$     xs = x + 1; %Add one to always be above zero.
    xs = x(find(x>0)); %Leave out all zero values.
end
%............................................................
%============================================================
%............................................................
n = length(xs);
%............................................................
xgms = exp(sum(log(xs))/n);
%............................................................
%============================================================
%............................................................
if all(x) %if all non-zero values.
    xgm = xgms;
else %if there were some zero values.
% $$$     xgm = xgms - 1;
    xgm = xgms;
end
%............................................................
%============================================================
%************************************************************
return

