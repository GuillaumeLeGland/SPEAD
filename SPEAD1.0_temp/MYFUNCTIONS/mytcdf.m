function p = tcdf(x,v);
%TCDF   Student's T cumulative distribution function (cdf).
%   P = TCDF(X,V) computes the cdf for Student's T distribution
%   with V degrees of freedom, at the values in X.
%
%   The size of P is the common size of X and V. A scalar input
%   functions as a constant matrix of the same size as the other input.

%   References:
%      [1] M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.7.
%      [2] L. Devroye, "Non-Uniform Random Variate Generation",
%      Springer-Verlag, 1986
%      [3] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, Section 10.3, pages 144-146.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 2.11 $  $Date: 2002/01/17 21:32:03 $

normcutoff = 1e7;
if nargin < 2,
    error('Requires two input arguments.');
end

[errorcode x v] = mydistchck(2,x,v);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
p=zeros(size(x));

% use special cases for some specific values of v
k = find(v==1);
    % See Devroye pages 29 and 450.
    % (This is also the Cauchy distribution)
if any(k)
    p(k) = .5 + atan(x(k))/pi;
end
k = find(v>=normcutoff);
if any(k)
    p(k) = normcdf(x(k));
end

% See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1
k = find(x ~= 0 & v ~= 1 & v > 0 & v < normcutoff);
if any(k),                            % first compute F(-|x|)
    xx = v(k) ./ (v(k) + x(k).^2);
    p(k) = betainc(xx, v(k)/2, 0.5)/2;
end

% Adjust for x>0.  Right now p<0.5, so this is numerically safe.
k = find(x > 0 & v ~= 1 & v > 0 & v < normcutoff);
if any(k), p(k) = 1 - p(k); end

p(x == 0 & v ~= 1 & v > 0) = 0.5;

% Return NaN for invalid inputs.
p(v <= 0 | isnan(x) | isnan(v)) = NaN;

