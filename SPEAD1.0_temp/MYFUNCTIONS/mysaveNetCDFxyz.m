function []=mymat2netCDF(X,VARNAME,rutaout)

%*********************************************
%Function MYMAT2NETCDF.m:
%
%This function saves a variable 3D(lat,long,time)=3D(180,360,12) 
%into a netCDF file.
%
%Use: mymat2netCDF(X,VARNAME)
%
%where:
%
%X: Variable 3D(180,360,12) we want to save.
%VARNAME: String with the 'name' of the variable.
%*********************************************    
    
%QUITO NAN:
Inan=find(isnan(X)==1);
X(Inan)=0;

%CREO UN FILE-netCDF:
filename=[VARNAME,'.nc'];
nc=netcdf(filename,'c');

%Define lengths:
nlatitude = [-89:1:90];
nlongitude = [-179:1:180];
ntime=[1:1:12];

%Define the dimension Time and Lat and Long of size 12 180 and 360 respectively.
nc('latitude')=180;
nc('longitude')=360;
nc('time')=12;

% coordinate variable longitudeitude
nc{'latitude'} = ncdouble('latitude');              % create a variable latitude of type double with 360 elements (dimension latitude).
nc{'latitude'}(:) = nlatitude;                       % store the octave variable latitude in the netcdf file
nc{'latitude'}.units = '[degrees]';                % define a string attribute of the variable latitude

% coordinate variable longitudeitude
nc{'longitude'} = ncdouble('longitude');              % create a variable longitude of type double with 360 elements (dimension longitude).
nc{'longitude'}(:) = nlongitude;                       % store the octave variable longitude in the netcdf file
nc{'longitude'}.units = '[degrees]';                % define a string attribute of the variable longitude

% coordinate variable time
nc{'time'} = ncdouble('time');;               % create a variable time of type double with 181 elements (dimension time). 
nc{'time'}(:) = ntime;                         % store the octave variable time in the netcdf file
nc{'time'}.units = '[months]';                % define a string attribute of the variable time

% variable DMS
nc{VARNAME} = ncdouble('latitude','longitude','time');        % create a variable DMS of type double of the size 360x180x12 (dimension lat,long,time).
nc{VARNAME}(:) = X;

%SAVE:
nc.history = 'netcdf file created by example_netcdf.m in octave'; %define a global string attribute
ncclose(nc)                              % close netcdf file and all changes are written to disk
disp(['example.nc file created. You might now inspect this file with the shell command "ncdump -h example.nc"']);

%******************************************
filename=[VARNAME,'.nc'];
file=[rutaout,filename];
if strcmp(rutaout,'junk')
    file=filename;
end
ncdump(file)
nc=netcdf(file,'nowrite');
VAR=nc{VARNAME}(:);
myglobalmaps(VAR,VARNAME,bone,0.02,'Monthly',1)
