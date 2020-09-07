function [Alowres] = mymeanbysegments1D(Ahigres,nbins)
%************************************************************************
%Use: [Alowres] = mymeanbysegments1D(Ahigres,nbins) 
%************************************************************************

%========================================================================
%........................................................................
nptos = length(Ahigres);
dx = nptos/nbins;
%........................................................................
Alowres = mean(reshape(Ahigres,dx,[]))';
%........................................................................
%========================================================================
%************************************************************************
return

% <http://www.mathworks.com/matlabcentral/newsreader/view_thread/319037>

% Create sample data.
A = randi(4, [9,1])

% Calculate mean of 3 rows, then the next 3 rows, etc.
jumping_mean = mean(reshape(A,3,[]))'

% Calculate means of rows 1-3, then of rows 2-4,
% then of rows 3-5, then of rows 4-6, etc.
sliding_mean = conv(A, [1;1;1]/3, 'valid')

