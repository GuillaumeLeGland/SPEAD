function m = geomean(x)
%GEOMEAN Geometric mean.
%   M = GEOMEAN(X) returns the geometric mean of the input.
%   When X is a vector with n elements, GEOMEAN(X) returns of the 
%   the n-th root of the product of the elements.
%   For a matrix input, GEOMEAN(X) returns a row vector containing
%   the geometric mean of each column of X.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 2.10 $  $Date: 2002/01/17 21:30:44 $

[r,n] = size(x);

% If the input is a row, make sure that N is the number of elements in X.
if r == 1, 
    r = n;
   x = x';
   n = 1; 
end

if any(any(x < 0))
    error('The data must all be non-negative numbers.')
end

m = zeros(1,n);
k = find(sum(x == 0) ==0);
m(k) = exp(sum(log(x(:,k))) ./ r);

