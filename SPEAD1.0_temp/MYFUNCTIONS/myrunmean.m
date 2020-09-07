function [varargout] = myrunmean(Xinput,runwin,varargin)
% $$$ function [Xsmooth,Xnoisy,Xstd]=myrunmean(Xinput,runwin,varargin)
% $$$ function [Xnoisy,Xsmooth,Xstd]=myrunmean(Xinput,runwin,varargin)
% $$$ function [Xsmooth,Xstd]=myrunmean(Xinput,runwin,varargin)
% $$$ function [Xsmooth,Xstd]=myrunmean(Xinput,runwin,wintype,keySignalFrequency)
% $$$ function [Xsmooth,Xstd]=myrunmean(Xinput,runwin,wintype)
% $$$ %function [Xsmooth,Xstd]=srunmean(Xinput,runwin,keySignalFrequency)
%-----------------------------------------------------------------------    
%Use:
%
% [Xsmooth,Xnoisy,Xstd]=myrunmean(Xinput,w,wintype,keySignalFrequency,keyWindowCentering)
%
% Xinput: vector o matriz que se quiere suavizar.
% runwin: vector [xwin,ywin,zwin] con el tamano de la mitad de la window deslizante (sin contar el pto central de ella; ej. si quiero window=3x3x3, runwin=[1,1,1]).
% 'wintype'= 'Espacial' / 'Espaciotemporal'.
% 'keySignalFrequency'= 'Periodica' / 'Aperiodica'.
% 'keyWindowCentering'= 'Centered'(default) / 'Backguard'.
%-----------------------------------------------------------------------
%(*) if 'keySignalFrequency'='Periodica' estoy diciendo que la senal es periodica y que
%puedo alargar los extremos con el final y el principio del vector (solo
%sirve para vectores). NO USAR!! (no se pq dije eso; pues para vectores parece que sale bien...)
%%-----------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SELECT HOW TO TREAT NAN VALUES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%........................................................................
% $$$ keyUseNANconv = 'not'; %PARA "globalPDR" PAPER USAR ESTE!!!
%........................................................................
keyUseNANconv = 'yes'; %Removes land pixels from the convolution.
%........................................................................

%%%%%%%%%%%%%%%%
%RUNMEAN WINDOW:
%%%%%%%%%%%%%%%%
%........................................................................
[msize,nsize,psize]=size(Xinput);
%........................................................................
% $$$ ArrayDim = ndims(Xinput); %NO USAR (does *not* give value of 1 for vectors)
ArrayDim = length(find(size(Xinput)-1) > 0); %OK (gives value of 1 for vectors)
%........................................................................
if ArrayDim==1
    xwinhalf = runwin(1);
    ywinhalf = [];
    zwinhalf = [];
    wintype = 'Espacial';
elseif ArrayDim==2
    xwinhalf = runwin(1);
    ywinhalf = runwin(2);
    zwinhalf = [];
    wintype = 'Espacial';
elseif ArrayDim==3
    xwinhalf = runwin(1); %Latitude.
    ywinhalf = runwin(2); %Longitude.
    zwinhalf = runwin(3); %Time.
    wintype = 'Espacial';
% $$$     wintype = 'Espaciotemporal'; %(gives some problems with ocean islands)
end
%........................................................................
xwin = (2*xwinhalf) + 1; %window size (it's always a odd number: 3,5,7,...etc)
ywin = (2*ywinhalf) + 1;
zwin = (2*zwinhalf) + 1;
%........................................................................
xyzwin = [xwin,ywin,zwin];
%........................................................................
winconv = ones(xyzwin); %window para la convolucion (matriz de 1's).
%........................................................................
% $$$ pcnt = 0; %minimo porcentage de puntos necesario para hacer el smooth.
pcnt = 1/3; %minimo porcentage de puntos necesario para hacer el smooth.
% $$$ % $$$ pcnt = 2/3; %minimo porcentage de puntos necesario para hacer el smooth.
% $$$ % $$$ pcnt = 3/4; %minimo porcentage de puntos necesario para hacer el smooth.
% $$$ % $$$ pcnt = 1.0; %minimo porcentage de puntos necesario para hacer el smooth.
%........................................................................
ptosmin = floor(pcnt*prod(xyzwin)); %que haya al menos un 50% de datos en la widow "xwin * ywin * zwin".
%........................................................................

%%%%%%%%%%%%%%%%%%%%%
%VARARGIN SELECTIONS:
%%%%%%%%%%%%%%%%%%%%%
numvarargin=length(varargin);
keyWindowCentering='Centered';
if numvarargin == 1
    wintype=varargin{1};
    keySignalFrequency='junk';
    keyWindowCentering='Centered';
elseif numvarargin == 2
    wintype=varargin{1};
    keySignalFrequency=varargin{2};
    keyWindowCentering='Centered';
elseif numvarargin == 3
    wintype=varargin{1};
    keySignalFrequency=varargin{2};
    keyWindowCentering=varargin{3};
end

%%%%%%%%%%%%%
%MATRIX SIZE:
%%%%%%%%%%%%%    
M = [msize,nsize,psize];
%................
% $$$ maxsize(M) %CREO QUE ASI NO FUNCIONA ENCONTRAR UN VECTOR, POR EJEMPLO SI [m,n,p]=[1,350,1]!!!!
%................
I = find(M==1);
%................
% $$$ if maxsize==1 %Xinput es un vector.
if length(I) >= 2 %Xinput es un vector.
    nlength = length(Xinput);
    X = Xinput(:);
    if numvarargin > 1
	if strcmp(keySignalFrequency,'Periodica')
	    %Si la senal es periodica le puedo añadir a los extremos el final
	    %y principio del vector, asi no cometo error en los extremos en
	    %el smoothing:
	    deltax=round(0.25*nlength);
	    X=[X(end-deltax:end);X;X(1:deltax)];
	    nlength=length(X);
	end
    end
    Xsmooth=ones(nlength,1)*nan;
    Xstd=ones(nlength,1)*nan;
    for i=1:nlength
	if i>xwin & i<nlength-xwin
	    xavei = mynanmean(X(i-xwin:i+xwin));
	    xstdi = mynanstd(X(i-xwin:i+xwin));
	elseif i<=xwin
	    xavei = mynanmean(X(1:i+xwin));
	    xstdi = mynanstd(X(1:i+xwin));
	elseif i >= (nlength-xwin)
	    xavei = mynanmean(X(i-xwin:nlength));
	    xstdi = mynanstd(X(i-xwin:nlength));
	end
	%.......................
	if strcmp(keyWindowCentering,'Centered')
	    Xsmooth(i) = xavei;
	    Xstd(i)    = xstdi;
	elseif strcmp(keyWindowCentering,'Backguard') %Si hago runmean de la semana anterior.
	    if (i+xwin) <= nlength
		Xsmooth(i+xwin) = xavei;
		Xstd(i+xwin)    = xstdi;
	    end
	end
	%.......................
    end
    if numvarargin > 1
	if strcmp(keySignalFrequency,'Periodica')
	    disp('periodic signal')
	    Xsmooth = Xsmooth(1+deltax+1:end-deltax); %quito los extremos que anadi.
	    Xstd = Xstd(1+deltax+1:end-deltax);
	end
    end
% $$$     [[1:n]',Xsmooth(:)]
% $$$     pause

    %%%%%%%%
    %OUTPUT:
    %%%%%%%%
    Xnoisy  = X(:)';
    Xsmooth = Xsmooth(:)'; %vector fila.
    Xstd    = Xstd(:)';

else %X es una matriz (2D o 3D)

    X = Xinput;

    if strcmp(keyUseNANconv,'not')
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%1. PONGO CEROS DONDE HAY -999, ZEROS and NAN's:
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%...............................
	Jneg = find(X < 0); 
	X(Jneg) = nan; %pongo NaN en las pos con -999.
	%...............................
	Izero = find(X < sqrt(eps) & X >= 0);
	X(Izero) = nan; %pongo NaN en las pos con casi-zero.
	%...............................
	Inan = find(isnan(X) == 1); 
	X(Inan) = 0; %pongo 0 en las pos con Nan.
	%...............................
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%3. DEFINO UNA MATRIZ CON UNOS EN DONDE HAY VALOR Y CEROS EN EL RESTO
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

    %%%%%%%%%%%%%%%%%
    %5. CONVOLUCION :   
    %%%%%%%%%%%%%%%%%
    if strcmp(wintype,'Espacial')
	for k=1:psize %psize(12meses)
	    %===========================================
	    %CONVOLUCION SOBRE LAS MATRICES CON VALORES:
	    %===========================================
	    if strcmp(keyUseNANconv,'not')
		%.............................................
		CX (:,:,k) = conv2(X(:,:,k),winconv(:,:,1),'same');             %E(X)=CX*(1/N)
		CXX(:,:,k) = conv2(X(:,:,k).*X(:,:,k),winconv(:,:,1),'same');   %E(Xexp2)=CXX*(1/N)
		%.............................................
	    elseif strcmp(keyUseNANconv,'yes')
		%.............................................
		CX (:,:,k) = mynanconv2(X(:,:,k),winconv(:,:,1),'same');             %E(X)=CX*(1/N)
		CXX(:,:,k) = mynanconv2(X(:,:,k).*X(:,:,k),winconv(:,:,1),'same');   %E(Xexp2)=CXX*(1/N)
		%.............................................
% $$$ 		CX (:,:,k) = mynanconv(X(:,:,k),winconv(:,:,1),'same');             %E(X)=CX*(1/N)
% $$$ 		CXX(:,:,k) = mynanconv(X(:,:,k).*X(:,:,k),winconv(:,:,1),'same');   %E(Xexp2)=CXX*(1/N)
		%.............................................
	    end
	    %========================================================================
	    %CONVOLUCION SOBRE LAS MATRICES CON 0 (donde Tierra) y 1 (donde valores):
	    %(es solo para contabilizar la cantidad de ptos usados en las convols)
	    %========================================================================
	    %.............................................
	    N(:,:,k) = conv2(X01(:,:,k),winconv(:,:,1),'same'); %n0 de ptos usados en la convolucion.
	    %.............................................
	end
    elseif strcmp(wintype,'Espaciotemporal')
	%....................
	%If size(180,360,12):
% $$$ 	twin=3; %si quiero promediar 3 meses.
	%....................
	%If size(180,360,365):
% $$$ 	twin=5; %si quiero promediar 5 dias.
% $$$ 	twin=7;
% $$$ 	twin=10;
	%....................

	tmpX01(:,:,1)           = X01(:,:,end); %anado 2 meses mas (Dec ano anterior y Jan ano proximo) para tener C.Frontera.
	tmpX01(:,:,1+psize+1)   = X01(:,:,1);
	tmpX01(:,:,1+1:1+psize) = X01;

	tmpX(:,:,1)           = X(:,:,end); %anado 2 meses mas (Dec ano anterior y Jan ano proximo) para tener C.Frontera.
	tmpX(:,:,1+psize+1)   = X(:,:,1);
	tmpX(:,:,1+1:1+psize) = X;
	
	%===========================================
	%CONVOLUCION SOBRE LAS MATRICES CON VALORES:
	%===========================================
	if strcmp(keyUseNANconv,'not')
	    %.............................................
	    CX  = convn(tmpX,      winconv,'same');   %E(X)=CX*(1/N)
	    CXX = convn(tmpX.*tmpX,winconv,'same');   %E(Xexp2)=CXX*(1/N)
	    %.............................................
	elseif strcmp(keyUseNANconv,'yes')
	    %.............................................
	    CX  = nanconvn(tmpX,      winconv,'same');   %E(X)=CX*(1/N)
	    CXX = nanconvn(tmpX.*tmpX,winconv,'same');   %E(Xexp2)=CXX*(1/N)
	    %.............................................
	end

	%========================================================================
	%CONVOLUCION SOBRE LAS MATRICES CON 0 (donde Tierra) y 1 (donde valores):
	%(es solo para contabilizar la cantidad de ptos usados en las convols)
	%========================================================================
	%.............................................
	N = convn(tmpX01,winconv,'same'); %n0 de ptos usados en la convolucion.
	%.............................................

	%==============================================
	%REMOVE THE PREVIOUSLY ADDED BORDER CONDITIONS:
	%==============================================
	N   = N  (:,:,2:psize+1); %quito las C.F.
	CX  = CX (:,:,2:psize+1);
	CXX = CXX(:,:,2:psize+1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %OBTENICION DE LAS MEDIAS MOBILES (SMOOTHING) Y  DE LAS VARIANZAS ASOCIADAS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-----------------------------
    %NOTA: VAR(X)=E(X^2)-(E(X))^2
    %-----------------------------
    %...............................
    In = find(N >= ptosmin); %pos(i,j,k) donde entraron mas de 33% de ptos en la window.
    %...............................
    EX  = ones(msize,nsize,psize)*nan;  %E(X).
    EX2 = ones(msize,nsize,psize)*nan;  %E(X^2).
    VX  = ones(msize,nsize,psize)*nan;  %Var(X).
    SX  = ones(msize,nsize,psize)*nan;  %Std(X).
    %...............................
    
    EX (In) = CX (In)./N(In); %E(X) = valores promdio en la window para el X
    EX2(In) = CXX(In)./N(In); %E(X^2)
    
    VX(In) = EX2(In) - EX(In).^2; %var(Xsmooth)
    VX(In) = fix(VX(In)*10^8)/10^8; %para evitar num negs debido a problemas de precision.
    SX(In) = sqrt(VX(In));

% $$$ Fc=N(In)./(N(In)-1);%Factor conversion Var con "N" grad lib, a Var con "N-1" grad lib.
% $$$ VX(In)=VX(In).*Fc;
% $$$ VY(In)=VY(In).*Fc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PONGO LOS NaN's DONDE ESTABAN:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %...............................
% $$$     X (Inan)=0;
% $$$     EX(Inan)=0;
% $$$     SX(Inan)=0;
    %...............................
    X (Inan)=nan; %PARA "globalPDR" PAPER USAR ESTOS!!!
    EX(Inan)=nan;
    SX(Inan)=nan;
    %...............................
% $$$     X (Inan)=-999;
% $$$     EX(Inan)=-999;
% $$$     SX(Inan)=-999;
    %...............................

    %%%%%%%%
    %OUTPUT:
    %%%%%%%%
    %===============================
    %SMOOTHED DATA POINTS AND NAN/ZERO VALUE DATA OTHERWISE:
    %...............................
    Xnoisy  = X; %PARA "globalPDR" PAPER USAR ESTOS!!!
    Xsmooth = EX;
    Xstd    = SX;
    %...............................
    %===============================
    %SMOOTHED DATA POINTS AND ORIGINAL NOISY DATA OTHERWISE:
    %...............................
% $$$     Xnoisy  = X;
% $$$     Xsmooth = Xnoisy; 
% $$$     Xstd    = SX; 
% $$$     %...............................
% $$$     Xsmooth(In) = EX(In); 
% $$$     %...............................
    %===============================
end

%%%%%%%%%%%%%%%%%%%%%
%DECREASE RESOLUTION:
%%%%%%%%%%%%%%%%%%%%%
%===================================================
%...................................................
% $$$ isub = [xwin:xwin:msize]; %lat
% $$$ jsub = [ywin:ywin:nsize]; %lon
% $$$ %%ksub = [zwin:zwin:psize]; %time
% $$$ ksub = [1:psize]; %time
%...................................................
isub = [2:2:msize]; %lat
jsub = [2:2:nsize]; %lon
ksub = [1:1:psize]; %time
%...................................................
[JSUB,ISUB,KSUB] = meshgrid(jsub,isub,ksub); %(cols,rows,planes)
%...................................................
H = sub2ind(size(X),ISUB,JSUB,KSUB);
%...................................................
XsmoothLowRes = Xsmooth(H);
XnoisyLowRes  = Xnoisy(H);
XstdLowRes    = Xstd(H);
%...................................................
%===================================================

%%%%%%%%%%%%%%
%FINAL OUTPUT:
%%%%%%%%%%%%%%
%...................................................
varargout{1} = Xsmooth;
varargout{2} = Xnoisy;
varargout{3} = Xstd;
%...................................................
varargout{4} = XsmoothLowRes;
varargout{5} = XnoisyLowRes;
varargout{6} = XstdLowRes;
%...................................................

%******************************************
return
%===================================================
%nanconv.m
%conv2nan.m
%...................................................
% <http://lasp.colorado.edu/cism/CISM_DX/code/CISM_DX-0.50/required_packages/octave-forge/extra/NaN/conv2nan.m>
% <http://sachinashanbhag.blogspot.com.es/2012/09/setting-up-random-number-generator-seed.html>
% <http://www.walkingrandomly.com/?p=2945>
%...................................................
%%rand('seed',24*03*1976) %random seed (old way) - DO *NOT* USE THIS WAY ANYMORE (see web post above)
rng(24*03*1976) %random seed (new way) - OKAY!
Xrandom = rand(5,7);
X = Xrandom;

X(3,4) = nan;
Y = ones(3,3); 
[Xconv] = conv2(X,Y,'same'); 

Inonan = ~isnan(X);
Jnonan = ~isnan(Y);

X(~Inonan) = 0;
Y(~Jnonan) = 0;

XN = real(Inonan);
YN = real(Jnonan);

Xones = ones(size(X));
Yones = ones(size(Y));

C = conv2(X,Y,'same');    % 2-dim convolution
N = conv2(XN,YN,'same');  % normalization term
F = conv2(Xones,Yones,'same'); % correction of normalization

Xnanconv = (C./N).*F
Nnannorm = N./F

XnanconvBis = C.*(F./N)
NnannormBis = 1./(F./N)

%===================================================
%nanconvn.m
%...................................................
rng(24*03*1976) %random seed (new way) - OKAY!
Xrandom = rand(5,7,3);
X = Xrandom;

X(3,4,2) = nan;
Y = ones(3,3,3); 
[Xconv] = convn(X,Y,'same'); 

Inonan = ~isnan(X);
Jnonan = ~isnan(Y);

X(~Inonan) = 0;
Y(~Jnonan) = 0;

XN = real(Inonan);
YN = real(Jnonan);

Xones = ones(size(X));
Yones = ones(size(Y));

C = convn(X,Y,'same');    % n-dim convolution
N = convn(XN,YN,'same');  % normalization term
F = convn(Xones,Yones,'same'); % correction of normalization

Xnanconv = (C./N).*F
Nnannorm = N./F

XnanconvBis = C.*(F./N)
NnannormBis = 1./(F./N)

%...................................................
%===================================================
