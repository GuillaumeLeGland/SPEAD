function [Apcnt_ari,Apcnt_geo] = myPercentageGeometric(A)

%========================================================================
%........................................................................
[msize,nsize,psize] = size(A); 
%........................................................................
ndim = ndims(A); 
Atot_ari = sum(A,ndim); 
Atot_geo = prod(A,ndim); 
Aave_geo = Atot_geo.^(1/psize);
%........................................................................
Apcnt_ari = A ./ repmat(Atot_ari,[1 1 psize]);
%%Apcnt_geo = A ./ repmat(Atot_geo,[1 1 psize]);
Apcnt_geo = A ./ repmat(Aave_geo,[1 1 psize]);
%........................................................................
sumApcnt_ari = sum(Apcnt_ari,3);
sumApcnt_geo = prod(Apcnt_geo,3);
%........................................................................
figure(200)
subplot(2,2,1)
imagesc(Atot_ari)
mycolorbar
subplot(2,2,2)
imagesc(Atot_geo)
mycolorbar
subplot(2,2,3)
imagesc(sumApcnt_ari)
mycolorbar
subplot(2,2,4)
imagesc(sumApcnt_geo)
mycolorbar
%........................................................................
figure(210)
counter = 0; 
for k = 1:psize
    counter =  counter + 1;
    Apcnt_arik = Apcnt_ari(:,:,k);
    Apcnt_geok = Apcnt_geo(:,:,k);
    subplot(2,3,counter)
    imagesc(Apcnt_arik)
    mycolorbar
    subplot(2,3,counter+3)
    imagesc(Apcnt_geok)
    mycolorbar
end
%........................................................................
%========================================================================
%************************************************************************
return

