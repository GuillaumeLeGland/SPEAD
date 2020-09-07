function [mymapjet]=mycolormapjet(varargin)

%*********************************************************************************
%Fucntion "MYCOLORMAPJET.m":
%---------------------------------------------------------------------------------
%This function gives you the "jet" colormap of the desired number of rows. 
%It takes the original [64 x 3] array, and transforms it into a [nbins x 3] array.
%Where "nbins" is the selected number of raws (eg. 12).
%---------------------------------------------------------------------------------
%
%Use: [mymapjet]=mycolormapjet(nbins)    
%
%
%*********************************************************************************

%%%%%%%%%%%%%%%%%%
%GET JET COLORMAP:
%%%%%%%%%%%%%%%%%%
mapjet=(jet);

%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE NUMBER OF BINS:
%%%%%%%%%%%%%%%%%%%%%%%
if nargin==1 %Manually.
    nbins=varargin{1};
else %Atomatically.
    [msize1,nsize1]=size(mapjet);
    nbins = msize1;
end

%%%%%%%%%%%%%%%
%ADD BLACK END:
%%%%%%%%%%%%%%%
%......................
% $$$ MAP1(1,:) = [0 0 0.15];
% $$$ MAP1(2,:) = [0 0 0.20];
% $$$ MAP1(3,:) = [0 0 0.25];
% $$$ MAP1(4,:) = [0 0 0.30];
% $$$ MAP1(5,:) = [0 0 0.35];
% $$$ MAP1(6,:) = [0 0 0.40];
% $$$ MAP1(7,:) = [0 0 0.45];
% $$$ MAP1(8,:) = [0 0 0.50];
% $$$ %......................
% $$$ MAP = [MAP1;mapjet];
%......................
MAP = [mapjet]; %NOT black end.
%......................

%%%%%%%%%%%%%
%INTERPOLATE:
%%%%%%%%%%%%%
[msize,nsize]=size(MAP);
x = [1:msize];
xmin=min(x);
xmax=max(x);
deltax=(xmax-xmin)/(nbins-1);
xi = [xmin:deltax:xmax];
[MAPI]=interp1(x,MAP,xi);
mymapjet = MAPI;


