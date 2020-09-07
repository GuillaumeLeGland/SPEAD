function [h,hcbar] = mynanimagesc(A)
%%function [h,hcbar] = mynanimagesc(A,cmap,nancolor)
% IMAGESC with NaNs assigning a specific color to NaNs
% <http://stackoverflow.com/questions/8481324/contrasting-color-for-nans-in-imagesc>

%# find minimum and maximum
amin=min(A(:));
amax=max(A(:));

%# size of colormap
cmap = jet; 
n = size(cmap,1);

%# color step
dmap=(amax-amin)/n;

%# standard imagesc
himg = imagesc(A);

%# add nan color to colormap
%%nancolor = [0.5 0.5 0.5]; %gray
%%nancolor = [1 1 1]; %white
nancolor = [0 0 0]; %black
%%nancolor = [0 1 1]; %cyan
colormap([nancolor; cmap]);

%# changing color limits
caxis([amin-dmap amax]);

%# place a colorbar
%%hcbar = colorbar;
hcbar = mycolorbar;

%# change Y limit for colorbar to avoid showing NaN color
ylim(hcbar,[amin amax])

if nargout > 0
    h = himg; %image handle.
end

%********************************************************
return
