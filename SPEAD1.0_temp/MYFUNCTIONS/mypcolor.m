function [hsurf] = pcolor(varargin)
%PCOLOR Pseudocolor (checkerboard) plot.
%   PCOLOR(C) is a pseudocolor or "checkerboard" plot of matrix C.
%   The values of the elements of C specify the color in each
%   cell of the plot. In the default shading mode, 'faceted',
%   each cell has a constant color and the last row and column of
%   C are not used. With shading('interp'), each cell has color
%   resulting from bilinear interpolation of the color at its 
%   four vertices and all elements of C are used. 
%   The smallest and largest elements of C are assigned the first and
%   last colors given in the color table; colors for the remainder of the 
%   elements in C are determined by table-lookup within the remainder of 
%   the color table.
%
%   PCOLOR(X,Y,C), where X and Y are vectors or matrices, makes a
%   pseudocolor plot on the grid defined by X and Y.  X and Y could 
%   define the grid for a "disk", for example.
%
%   PCOLOR(AX,..) plots into AX instead of GCA.
%
%   H = PCOLOR(...) returns a handle to a SURFACE object.
%
%   PCOLOR is really a SURF with its view set to directly above.
%
%   See also CAXIS, SURF, MESH, IMAGE, SHADING.

%-------------------------------
%   Additional details:
%
%
%   PCOLOR sets the View property of the SURFACE object to directly 
%   overhead.
%
%   If the NextPlot axis property is REPLACE (HOLD is off), PCOLOR resets 
%   all axis properties, except Position, to their default values
%   and deletes all axis children (line, patch, surf, image, and 
%   text objects).  View is set to [0 90].

%   Copyright 1984-2006 The MathWorks, Inc. 
%   $Revision: 5.9.4.6 $  $Date: 2011/07/25 03:49:30 $

%   J.N. Little 1-5-92

% Parse possible Axes input
[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,3,nargs,'struct'))

% do error checking before calling newplot. This argument checking should
% match the surface(x,y,z) or surface(z) argument checking.
if nargs == 2
  error(message('MATLAB:pcolor:InvalidNumberOfInputs'))
end
if isvector(args{end})
  error(message('MATLAB:pcolor:NonMatrixColorInput'));
end
if nargs == 3 && LdimMismatch(args{1:3})
  error(message('MATLAB:pcolor:InputSizeMismatch'));
end
for k = 1:nargs
  if ~isreal(args{k})
    error(message('MATLAB:pcolor:NonRealInputs'));
  end
end

cax = newplot(cax);
hold_state = ishold(cax);
%========================================================================
%........................................................................
if nargs == 1 %Padding works OKAY.
    Array = args{1};
    padval = mean(Array(:)); %use mean value for padding extra row and column.
    Arraypad = padarray(Array,[1 1],padval,'post'); %Adding one extra row and column just for plotting like imagesc.
    hsurf = surface(zeros(size(Arraypad)),Arraypad,'parent',cax);
    [m,n] = size(Arraypad);
    lims = [ 1 n 1 m];
    set(gca,'Xlim',[1 n],'Xtick',[0:n-1]+0.5,'XtickLabel',[0:n-1]) %Centering labels.
    set(gca,'Ylim',[1 m],'Ytick',[0:m-1]+0.5,'YtickLabel',[0:m-1])
elseif nargs == 3 %Padding does NOT WORK HERE!!! 
    [X,Y,Array] = deal(args{1:3});
%%    padval = mean(Array(:)); %use mean value for padding extra row and column.
    padval = nan;
    Arraypad = padarray(Array,[1 1],padval,'post'); %Adding one extra row and column just for plotting like imagesc.
% $$$     Xpad = padarray(X,[1 1],'replicate','post'); %Adding one extra row and column just for plotting like imagesc.
% $$$     Ypad = padarray(Y,[1 1],'replicate','post'); %Adding one extra row and column just for plotting like imagesc.
    Xpad = padarray(X,[1 1],nan,'post'); %Adding one extra row and column just for plotting like imagesc.
    Ypad = padarray(Y,[1 1],nan,'post'); %Adding one extra row and column just for plotting like imagesc.
% $$$     Xpad,Ypad,Arraypad
% $$$     pause 
    hsurf = surface(Xpad,Ypad,zeros(size(Arraypad)),Arraypad,'parent',cax);
    lims = [min(min(X)) max(max(X)) min(min(Y)) max(max(Y))];
end
%........................................................................
%========================================================================
%........................................................................
keyMapOrientation = 'PCOLOR';
%%keyMapOrientation = 'IMAGESC';
%........................................................................
if strcmp(keyMapOrientation,'PCOLOR')
    set(gca,'ydir','normal') %For pcolor orientation.
elseif strcmp(keyMapOrientation,'IMAGESC')
    set(gca,'ydir','reverse') %For imagesc orientation.
end
%........................................................................
%========================================================================
%***************************************************
return

if ~hold_state
    set(cax,'View',[0 90]);
    set(cax,'Box','on');
    if lims(2) <= lims(1)
        lims(2) = lims(1)+1;
    end
    if lims(4) <= lims(3)
        lims(4) = lims(3)+1;
    end
    axis(cax,lims);
end
if nargout == 1
    h = hsurf;
end

function ok = LdimMismatch(x,y,z)
[xm,xn] = size(x);
[ym,yn] = size(y);
[zm,zn] = size(z);
ok = (xm == 1 && xn ~= zn) || ...
     (xn == 1 && xm ~= zn) || ...
     (xm ~= 1 && xn ~= 1 && (xm ~= zm || xn ~= zn)) || ...
     (ym == 1 && yn ~= zm) || ...
     (yn == 1 && ym ~= zm) || ...
     (ym ~= 1 && yn ~= 1 && (ym ~= zm || yn ~= zn));



