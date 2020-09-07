function []=mymonthlymeansbyoceanregionfigure(varargin)

numvar=length(varargin);
[ZONAS,TERRA,MAPterra]=myoceanregions;

if numvar==1
    [mX]=mymonthlymeansbyoceanregion(X,ZONAS);
    mymonthlymeansbyoceanregionplots(VARNAME,TERRA,MAPterra,1,mX);
elseif numvar==2
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
mymonthlymeansbyoceanregionplots(VARNAME,TERRA,MAPterra,1,mX);
elseif numvar==3
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
mymonthlymeansbyoceanregionplots(VARNAME,TERRA,MAPterra,1,mX);
elseif numvar==4
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
mymonthlymeansbyoceanregionplots(VARNAME,TERRA,MAPterra,1,mX);
elseif numvar==5
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
mymonthlymeansbyoceanregionplots(VARNAME,TERRA,MAPterra,1,mX);
elseif numvar==6
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
mymonthlymeansbyoceanregionplots(VARNAME,TERRA,MAPterra,1,mX);
elseif numvar==7
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
elseif numvar==8
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
end

[mX]=mymonthlymeansbyoceanregion(X,ZONAS);
mymonthlymeansbyoceanregionplots(VARNAME,TERRA,MAPterra,1,mX);
