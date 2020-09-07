function [EX,SX] = myrunaveconv2D(Xinput,xwin,ywin)

%========================================================================
%........................................................................
X = Xinput;
xyzwin = [xwin,ywin];
winconv = ones(xyzwin); %window para la convolucion (matriz de 1's).
pcnt = 1/3; %minimo porcentage de puntos necesario para hacer el smooth.
ptosmin = floor(pcnt*prod(xyzwin)); %que haya al menos un 33% de datos en la widow "xwin * ywin".
[msize,nsize,psize]=size(Xinput);
%........................................................................
%========================================================================
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SELECT HOW TO TREAT NAN VALUES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
% $$$ keyUseNANconv = 'not'; %PARA "globalPDR" PAPER USAR ESTE!!!
%........................................................................
keyUseNANconv = 'yes'; %Removes land pixels from the convolution.
%........................................................................
%========================================================================
if strcmp(keyUseNANconv,'not')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PONGO CEROS DONDE HAY -999, ZEROS and NAN's:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Jneg = find(X < 0); 
    X(Jneg) = nan; %pongo NaN en las pos con -999.

    Izero = find(X < sqrt(eps) & X >= 0);
    X(Izero) = nan; %pongo NaN en las pos con casi-zero.

    Inan = find(isnan(X) == 1); 
    X(Inan) = 0; %pongo 0 en las pos con Nan.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DEFINO UNA MATRIZ CON UNOS EN DONDE HAY VALOR Y CEROS EN EL RESTO
    %(TIERRA Y ZONAS SIN DATOS) PARA OBTENER EL NUMERO DE GRADOS DE LIBERTAD
    %(O SEA, CANTIDAD DE PTOS USADOS EN EL CALCULO DE LOS ESTADISTICOS):
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X01 = X;
    Ipos = find(X > 0);
    X01(Ipos) = 1;

elseif strcmp(keyUseNANconv,'yes')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %COUNT NUMBER OF DATA POINTS WITHOUT TAKING NANs:
    %[NOTE: Use with "mynanconv.m"]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X01 = zeros(size(X));
    Inan = find(isnan(X) == 1); %NAN data.
    Inonan = find(isfinite(X) == 1); %non-nan data (zero is a valid value).
    X01(Inonan) = 1; %Array with ones where valid data and zero otherwise.

    end
    
for k=1:psize %psize(12meses)

    %CONVOLUCION SOBRE LAS MATRICES CON VALORES:
    if strcmp(keyUseNANconv,'not')
	CX (:,:,k) = conv2(X(:,:,k),winconv(:,:,1),'same');             %E(X)=CX*(1/N)
	CXX(:,:,k) = conv2(X(:,:,k).*X(:,:,k),winconv(:,:,1),'same');   %E(Xexp2)=CXX*(1/N)
    elseif strcmp(keyUseNANconv,'yes')
	CX (:,:,k) = mynanconv2(X(:,:,k),winconv(:,:,1),'same');             %E(X)=CX*(1/N)
	CXX(:,:,k) = mynanconv2(X(:,:,k).*X(:,:,k),winconv(:,:,1),'same');   %E(Xexp2)=CXX*(1/N)
    end
    %========================================================================
    %CONVOLUCION SOBRE LAS MATRICES CON 0 (donde Tierra) y 1 (donde valores):
    %(es solo para contabilizar la cantidad de ptos usados en las convols)
    %========================================================================
    N(:,:,k) = conv2(X01(:,:,k),winconv(:,:,1),'same'); %n0 de ptos usados en la convolucion.
end
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENICION DE LAS MEDIAS MOBILES (SMOOTHING) Y  DE LAS VARIANZAS ASOCIADAS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------
%NOTA: VAR(X)=E(X^2)-(E(X))^2
%-----------------------------
%========================================================================
%........................................................................
In = find(N >= ptosmin); %pos(i,j,k) donde entraron mas de 33% de ptos en la window.
%........................................................................
EX  = ones(msize,nsize,psize)*nan;  %E(X).
EX2 = ones(msize,nsize,psize)*nan;  %E(X^2).
VX  = ones(msize,nsize,psize)*nan;  %Var(X).
SX  = ones(msize,nsize,psize)*nan;  %Std(X).
%........................................................................
EX (In) = CX (In)./N(In); %E(X) = valores promdio en la window para el X
EX2(In) = CXX(In)./N(In); %E(X^2)
%........................................................................
VX(In) = EX2(In) - EX(In).^2; %var(Xsmooth)
VX(In) = fix(VX(In)*10^8)/10^8; %para evitar num negs debido a problemas de precision.
SX(In) = sqrt(VX(In));
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PONGO LOS NaN's DONDE ESTABAN:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
X (Inan)=nan;
EX(Inan)=nan;
SX(Inan)=nan;
%........................................................................
%========================================================================

return
