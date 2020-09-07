function []=mysaveNetCDF(VAR,VARNAME,pathout,filenameout,varunits,varargin)

%Use: mysaveNetCDF(VAR,VARNAME,pathout,filenameout,varunits,varargin)
%
%<http://ecco2.jpl.nasa.gov/data1/matlab/toolbox3/netcdf/ncutility/ncexample.m>
% ncexample.m -- "NetCDF Toolbox for Matlab-5" example.
%  ncexample (no argument) is a short example that lists
%   itself, builds a simple NetCDF file, then displays
%   its variables.
 
% Copyright (C) 1997 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 12-Jun-1997 16:23:04.

type(mfilename)

help(mfilename)

%DIMENSIONS NAMES AND UNITS:
%.....................
%IF varargin==0:
xtag='junk';
xunits='junk';

ytag='junk';
yunits='junk';

ztag='junk';
zunits='junk';

wtag='junk';
wunits='junk';
%.....................
nvarargin=length(varargin);    
if nvarargin==2
    xtag=varargin{1}
    xunits=varargin{2}
elseif nvarargin==4
    xtag=varargin{1}
    xunits=varargin{2}
    
    ytag=varargin{3}
    yunits=varargin{4}
elseif nvarargin==6
    xtag=varargin{1}
    xunits=varargin{2}

    ytag=varargin{3}
    yunits=varargin{4}

    ztag=varargin{5}
    zunits=varargin{6}
elseif nvarargin==8
    xtag=varargin{1}
    xunits=varargin{2}

    ytag=varargin{3}
    yunits=varargin{4}

    ztag=varargin{5}
    zunits=varargin{6}

    wtag=varargin{7}
    wunits=varargin{8}
end
%.....................

%ARRAY SIZE:
[m,n,o,p]=size(VAR);
VARsize=[m,n,o,p]

%QUITO NAN:
Inan=find(isnan(VAR)==1);
VAR(Inan)=0;

%DEFINE AXIS NAMES:
axislabelX=['X.',xtag]; %ojo, NO dejar espacios entre X y xtab ('Xxtag' SI; 'X xtag' NO!!)
axislabelY=['Y.',ytag];
axislabelZ=['Z.',ztag];
axislabelW=['W.',wtag];

% ---------------------------- DEFINE THE FILE --------------------------- %
%No NetCDF warnings:
ncquiet                                              

%Create NetCDF file:
if strcmp(pathout,'base')
    fileout=filenameout;
else
    fileout=[pathout,filenameout];
end
fileout
%????????????????
% $$$ fileout='borrar.nc';
% $$$ fileout
%????????????????
nc=netcdf(fileout,'clobber');

%Global attributes:
nc.description = 'Creating NetCDF file';
nc.author = 'Dr. Sergio M. Vallina';
%nc.date = '31 Jan 2008';
%nc.date = '2 Apr 2008';
nc.date = date;

%Define dimensions:
nc(axislabelX)=m;
nc(axislabelY)=n;
nc(axislabelZ)=o;
nc(axislabelW)=p;

%Define variables:
nc{axislabelX} = axislabelX;
nc{axislabelY} = axislabelY;
nc{axislabelZ} = axislabelZ;
nc{axislabelW} = axislabelW;
nc{VARNAME} = {axislabelX,axislabelY,axislabelZ,axislabelW};

%Attributes.
nc{axislabelX}.units = xunits;
nc{axislabelY}.units = yunits;
nc{axislabelZ}.units = zunits;
nc{axislabelW}.units = wunits;
nc{VARNAME}.units = varunits;

% ---------------------------- STORE THE DATA ---------------------------- %

%Define lengths:
nx=[1:m];
ny=[1:n];
nz=[1:o];
nw=[1:p];

%Put all the data.
nc{axislabelX}(:) = nx;
nc{axislabelY}(:) = ny;
nc{axislabelZ}(:) = nz;
nc{axislabelW}(:) = nw;
nc{VARNAME}(:) = VAR;

nc = close(nc);                                      % Close the file.

% ---------------------------- RECALL THE DATA --------------------------- %
% $$$ %................................
% $$$ file=[pathout,filenameout];
% $$$ if strcmp(pathout,'junk')
% $$$     file=filenameout;
% $$$ end
% $$$ ncdump(file)
% $$$ %Open NetCDF file.
% $$$ nc=netcdf(file,'nowrite');
% $$$ %GET VARIABLE:
% $$$ VARout=nc{VARNAME}(:);
% $$$ %PLOT VARIABLE:
% $$$ myglobalmaps(VARout,VARNAME,bone,0.02,'Monthly',1)
% $$$ %................................
% $$$ nc = netcdf('myncexample.nc', 'nowrite');              % Open NetCDF file.
% $$$ description = nc.description(:)                      % Global attribute.
% $$$ variables = var(nc);                                 % Get variable data.
% $$$ for i = 1:length(variables)
% $$$    disp([name(variables{i}) ' =']), disp(' ')
% $$$    disp(variables{i}(:))
% $$$ end
% $$$ nc = close(nc);                                      % Close the file.
% $$$ %................................
% --------------------------------- DONE --------------------------------- %
