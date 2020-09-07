function [GLOBE,Landgcm] = myLandMaskMITgcm()

topo = load('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MODELLING/MITesm/MITesm3D/DATA-OFFLINE/INITIAL-CONDITIONS/ASCII/gcmLAND.dat');

mlat = 180;
nlon = 360;

depthz0 = 0; %[m]

Landgcm = find(topo >= depthz0);

% $$$ GLOBE = zeros(mlat,nlon);
% $$$ GLOBE(Landgcm) = 1.0;

% $$$ GLOBE = ones(mlat,nlon);
% $$$ GLOBE(Landgcm) = zeros;

GLOBE = ones(mlat,nlon);
GLOBE(Landgcm) = nan;

