function [M] = mysequentialmatrix(nx,ny)

%**********************************************
%Use [Matrix] = mysequentialmatrix(ncols,nrows)
%**********************************************

M = reshape([1:nx*ny],[nx ny])';

return

