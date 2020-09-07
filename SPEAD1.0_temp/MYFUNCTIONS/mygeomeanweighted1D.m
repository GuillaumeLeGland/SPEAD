function [gsm,gsd] = mygeomeanweighted1D(x,w)

%************************************************************
% Use: [gsm,gsd] = mygeomeanweighted1D(x,w)
%************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK FOR NEGATIVE AND ZERO DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%============================================================
%............................................................
if any(x < 0)
    error('The data must all be non-negative numbers.')
end
%............................................................
if all(x) %if all non-zero values.
    xs = x; 
    ws = w;
else %if there are some zero values.
    Inonzero = find(x > 0);
    xs = x(Inonzero); %Leave out all zero values.
    ws = w(Inonzero);
end
%............................................................
%============================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GEOMETRIC SAMPLE MEAN AND STANDARD DEVIATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------
% <https://www.researchgate.net/post/How_can_I_calculate_the_value_of_the_geometric_standard_deviation_taking_into_account_weight>
%------------------------------------------------------------
%............................................................
gsm = exp(sum(ws.*log(xs))/sum(ws)); %Geometric sampling mean.
gsv = sum(ws .* (log(xs)-log(gsm)).^2)/sum(ws); %Geometric sampling variance.
gsd = exp(sqrt(gsv)); %Geometric sampling deviation.
%............................................................
%============================================================
%************************************************************
return

