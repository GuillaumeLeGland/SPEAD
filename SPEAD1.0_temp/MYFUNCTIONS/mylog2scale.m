function [Xstar]=mylog2scale(X)

%.........................
xmin = min(X(:));
xmax = max(X(:));
%.........................
Log2xmin = log2(xmin);
Log2xmax = log2(xmax);
%.........................
nbins = 256;
%.........................
alfa =  Log2xmin;
beta = (Log2xmax - Log2xmin) / nbins;
%.........................
Xlog2 = (log2(X) - alfa) / beta; %[0 - 256]
%.........................
nXlog2 = Xlog2 / nbins; %[0 - 1]
nX     = X     / xmax;
%.........................
if ndims(X) == 2 %Arrays 2D
    Xstar(:,:,1) = X;
    Xstar(:,:,2) = nX;
    Xstar(:,:,3) = nXlog2;
elseif ndims(X) == 1 %Does NOT recognize vectors (so it's useless).
    Xstar(1,:) = X;
    Xstar(2,:) = nX;
    Xstar(3,:) = nXlog2;
end
%.........................
Xstar = squeeze(Xstar);
%.........................
