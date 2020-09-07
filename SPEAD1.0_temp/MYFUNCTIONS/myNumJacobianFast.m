function [Jf]=myNumJacobianFast(f,X) 

%**********************************
%FUNCTION: myNumJacobianFast.m
%Approximately calculate Jacobian matrix.
%Source: <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>
%**********************************
m  = length(X);
Jf = zeros(m,m);
dx = 0.001;
for j = 1:m
    %......................
    X1    = X;
    X2    = X;
    %......................
    X2(j) = X(j) + dx;
    %......................
% $$$     f1 = f(X1); %This way does NOT work in Matlab6p5.
% $$$     f2 = f(X2);
    %......................
    f1 = feval(f,X1); %This way allways works ok.
    f2 = feval(f,X2);
    %......................
    Jf(:,j) = (f2 - f1)/dx;  % partial derivatives in j-th row.
    %......................
end
%%Jf = Jf'; %Do NOT use.
