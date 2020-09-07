function [Xnanconv] = mynanconv2(X,Y,shape)
% <http://lasp.colorado.edu/cism/CISM_DX/code/CISM_DX-0.50/required_packages/octave-forge/extra/NaN/conv2nan.m>

Inonan = ~isnan(X);
Jnonan = ~isnan(Y);

X(~Inonan) = 0;
Y(~Jnonan) = 0;

XN = real(Inonan);
YN = real(Jnonan);

Xones = ones(size(X));
Yones = ones(size(Y));

C = conv2(X,Y,shape);    % 2-dim convolution
N = conv2(XN,YN,shape);  % normalization term
F = conv2(Xones,Yones,shape); % correction of normalization

Xnanconv = C.*(F./N);
Nnannorm = 1./(F./N);



