function [DistancePoint1vsPoint2]=myDistanceBetween2points(point1,point2)
% $$$ function [DistancePoint1vsPoint2,RelativeDistancePoint1vsPoint2]=myDistanceBetween2points(point1,point2)

%***********************************
%Programme myDistanceBetween2points.m:
%This script calculates the geometric distance between two points in a
%n-dimensional space using the formula:
%
% point1(x1,y1,z1)
% point2(x2,y2,z2)
%
% distance = sqrt( (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 )
%
%   
%Use:
%
% [DistancePoint1vsPoint2]=myDistanceBetween2points(point1,point2)
%
%***********************************
%.......................
nsize1 = length(point1);
nsize2 = length(point2);
%.......................
if nsize1~=nsize2
    display('error! point1 and point2 have to be of same length!')
    pause
else
    n=nsize1; %space dimension.
end
%.......................

%============
%Remove NaNs:
%============
%........................
Inan1 = isnan(point1)==0; %Logical vector with ONES in the elementes with numbers, and ZERO in the elemens with NaNs.
Inan2 = isnan(point2)==0;
%........................
Inan = Inan1.*Inan2; %ZEROS where either point1 or point2 had NaNs.
%........................
Inan = logical(Inan); %convert to "logical" vector again, so that it can be used for indexing.
%........................
point1 = point1(Inan); %Vector without NaNs.
point2 = point2(Inan); %Vector without NaNs.
%........................

%=============================
%Distance for each coordenate:
%=============================
%........................
msize = length(point1); %Space dimensions.
%........................
for j=1:msize
    %............
    DistAbs(j)=(point2(j)-point1(j))^2; %Absolute distance (vector modulus) %USAR ESTE!!!!
    %............
% $$$     DistRel(j)=((point2(j)-point1(j))/point1(j))^2; %Relative distance.
    %............
end

%==============================
%Distance between the 2 points:
%==============================
DistancePoint1vsPoint2Abs = sqrt(sum(DistAbs));
% $$$ DistancePoint1vsPoint2Rel = sqrt(sum(DistRel));

%=======
%OUTPUT:
%=======
DistancePoint1vsPoint2=DistancePoint1vsPoint2Abs;
% $$$ DistancePoint1vsPoint2=DistancePoint1vsPoint2Rel;
