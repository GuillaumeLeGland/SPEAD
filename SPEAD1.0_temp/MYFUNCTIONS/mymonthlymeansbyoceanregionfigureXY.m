function []=mymonthlymeansbyoceanregionfigureXY(VARNAMEx,VARNAMEy,X,Y,xscalemax,yscalemax,dx,dy,fignum)

%**************************************
%Use: []=mymonthlymeansbyoceanregionfigureXY('VARNAMEx','VARNAMEy',X,Y,xscalemax,yscalemax,dx,dy)
%**************************************
    
[ZONAS,TERRA,MAPterra]=myoceanregions;

[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mY]=mymonthlymeansbyoceanregion(Y,ZONAS);
mymonthlymeansbyoceanregionplotsXY(VARNAMEx,VARNAMEy,TERRA,MAPterra,fignum,mX,mY,xscalemax,yscalemax,dx,dy);
