function [VAR,ncinfo]=myreadNetCDF(filepath,filename,VARNAME)

%*********************************************
%Function MYREADNETCDF.m:
%
%This function reads a variable 3D(lat,long,time)=3D(180,360,12) 
%from a netCDF file.
%
%Use: [VAR]=myreadNetCDF(filepath,filename,VARNAME)
%
%where:
%
%VAR: variable 3D(180,360,12) we want to read.
%filepath: string with the '/path/' were is the netCDF file (type 'base' for current directory)
%filename: string with the 'filename.nc' of the netCDF file.
%VARNAME: string with the 'name' of the variable.
%*********************************************    
%(see also # <http://codes.guillaumemaze.org/tips/howtomakecleannetcdffilefrommatlab>)

%%%%%%%%%%%
%LOAD FILE:
%%%%%%%%%%%
filenetcdf = [filepath,filename];
if strcmp(filepath,'base')
    filenetcdf = filename;
    %%ncdump(filenetcdf)
end
display(filenetcdf)
nc = netcdf(filenetcdf,'nowrite');
%%dim(nc) 
%%coord(nc) 
VAR = nc{VARNAME}(:); 
ncinfo = ncvar(VARNAME,nc); 

