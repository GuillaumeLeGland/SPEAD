function [Ex,Ey,Exy,Vx,Vy,N]=myspatialcorrmaps_estadisticos(X,Y,h,L,mode) 
%'mode' debe ser un string definiendo si quiero correlacion mensual (12
%mapas de correlacion, uno para cada mes) o si quiero una correlacion
%anual (un solo mapa de correlacion integrando todos los puntos del ano).
    
A=X; %DMS
B=Y; %CCN
[m,n,p]=size(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. QUITO LOS NaN's Y PONGO ZEROS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A0=A; %defino nueva variable A0 sobre la que pondre ceros en la Tierra y no.data
B0=B;
Inan=find(isnan(X)==1);
Jnan=find(isnan(Y)==1);
A0(Inan)=0;
B0(Jnan)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. PONGO CEROS DONDE HAY NaN's O LAND:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VARX:
Ia=find(isnan(A)==1); %posic. con Nan.
Ja=find(A<0); %pos con -999
A0=A; %defino nueva variable A0 sobre la que pondre ceros en la Tierra y no.data
if length(Ia)>0
    A0(Ia)=0;
end
if length(Ja)>0
    A0(Ja)=0;
end
%VARY:
Ib=find(isnan(B)==1); %posic. con Nan.
Jb=find(B<0); %pos con -999
B0=B;
if length(Ib)>0
    B0(Ib)=0;
end
if length(Jb)>0
    B0(Jb)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3. PONGO CEROS EN LAS POS(i,j) DONDE A0 y/o B0 NO TIENEN VALORES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KA=find(A0==0);
KB=find(B0==0);
K=[KA;KB];

A0(K)=0;
B0(K)=0;
    
FA0=find(A0==0);
FB0=find(B0==0);

if FA0~=FB0, display('cuidado!! no esta bien hecho lo de poner ceros donde A y/o B no tengan valores')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. DEFINO UNA MATRIZ CON UNOS EN DONDE HAY VALOR Y CEROS EN EL RESTO
%(TIERRA Y ZONAS SIN DATOS) PARA OBTENER EL NUMERO DE GRADOS DE LIBERTAD
%(O SEA, CANTIDAD DE PTOS USADOS EN EL CALCULO DE LOS ESTADISTICOS):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A01=A0;
B01=B0;
G=find(A0>0 | B0>0);
A01(G)=1;
B01(G)=1;

%%%%%%%%%%%%%%%%%
%5. CONVOLUCION :   
%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------
%NOTA: La convolucion me permite obtener las medias y varianzas 
%con las que calculare el Coef. de correlacion de PEARSON:
% r=(E(xy)-E(x)E(y))/sqrt(Var(x)*Var(y)) 
%---------------------------------------------------------------
window=ones(h); %window para la convolucion (matriz de 1's)

for k=1:p %p=12meses
    %===========================================================================
    %1. CONVOLUCION SOBRE LAS MATRICES CON 0 (donde Tierra) y 1 (donde valores):
    %(es solo para contabilizar la cantidad de ptos usados en las convols)
    %============================================================================
    CA01(:,:,k)=conv2(A01(:,:,k),window,'same');
    CB01(:,:,k)=conv2(B01(:,:,k),window,'same');
    if CA01(:,:,k)~=CB01(:,:,k)
        display('Ojo!, el n0 de ptos usados en la convolucion no es el mismo para A y B')
            pause
    else 
        N(:,:,k)=CA01(:,:,k); %n0 de ptos usados en la convolucion.
    end
    
    %==============================================
    %2. CONVOLUCION SOBRE LAS MATRICES CON VALORES:
    %==============================================
    CA0(:,:,k)=conv2(A0(:,:,k),window,'same');                %E(A)=CA0*(1/N)
    CB0(:,:,k)=conv2(B0(:,:,k),window,'same');                %E(B)=CB0*(1/N)
    CA0A0(:,:,k)=conv2(A0(:,:,k).*A0(:,:,k),window,'same');   %E(Aexp2)=CA0A0*(1/N)
    CB0B0(:,:,k)=conv2(B0(:,:,k).*B0(:,:,k),window,'same');   %E(Bexp2)=CB0B0*(1/N)
    CA0B0(:,:,k)=conv2(A0(:,:,k).*B0(:,:,k),window,'same');   %E(AB)=CA0B0*(1/N)
end

%=========================================================
%3. COJO SOLO LA CONVOLUCION DONDE SE USO AL MENOS UN PTO (O SEA, PONGO
%CEROS DONDE NO SE USO NI SIQUIERA UN VALOR VALIDO):
%=========================================================
H=find(N==0);
CA0(H)=0;
CB0(H)=0;
CA0A0(H)=0;
CB0B0(H)=0;
CA0B0(H)=0;

if p>1 %X e Y son matrices 3D

    if strcmp(lower(mode(1:5)),'anual') %CORRELACION ANUAL

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %4. CALCULO LAS ESPERANZAS (MEDIA Y VARIANZA) ANUALES (SON PARA Rpearson):
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %===============================================================
        %1. SUMO LAS CONVOLUCIONES Y EL No DE PTOS USADO A LO LARGO DEL 
        %EJE TEMPORAL (meses):
        %===============================================================
        sumCA0=sum(CA0,3); 
        sumCB0=sum(CB0,3);
        sumCA0A0=sum(CA0A0,3);
        sumCB0B0=sum(CB0B0,3);
        sumCA0B0=sum(CA0B0,3);
        sumN=sum(N,3);
        sumCA01=sum(CA01,3);
        sumCB01=sum(CB01,3);
        
        %======================
        %2. OBTENGO LAS MEDIAS:
        %======================
        Am=ones(m,n)*nan;
        Bm=ones(m,n)*nan;
        A2m=ones(m,n)*nan;
        B2m=ones(m,n)*nan;
        ABm=ones(m,n)*nan;
        
        %ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
        %COJO SOLO LAS POS(i,j) DONDE SE TIENEN UN MINIMO DE VALORES
        %PREDEFINIDO, PARA EL CALCULO DE LA CORRELACION ANUAL:
        %ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
        if h==3
            I=find(sumN>=0.5*h*h*p);
            J=find(sumN< 0.5*h*h*p);
        elseif h>3
            I=find(sumN>=100);
            J=find(sumN< 100);
        end

        Am(I) = sumCA0(I)./sumN(I);       %E(A) %Am(55:65,200:210)
        Bm(I) = sumCB0(I)./sumN(I);       %E(B)
        A2m(I)= sumCA0A0(I)./sumN(I);     %E(A2)
        B2m(I)= sumCB0B0(I)./sumN(I);     %E(B2)
        ABm(I)= sumCA0B0(I)./sumN(I);     %E(AB)

        %=============================================
        %3. OBTENGO LAS VARIANZAS: VAR(X)=E(X2)-E2(X)
        %=============================================
        Av=A2m-(Am.^2);    %Var(A)
        Bv=B2m-(Bm.^2);    %Var(B)

% $$$         Fc=sumN(I)./(sumN(I)-1);%Fact.conv.pasar d Var.con "N" grad.lib, a Var.con "N-1" grad.lib.
% $$$         Av(I)=Av(I).*Fc;
% $$$         Bv(I)=Bv(I).*Fc;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PONGO LA MASCARA TERRESTRE:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Am(L)=nan;
        Bm(L)=nan;
        ABm(L)=nan;
        Av(L)=nan;
        Bv(L)=nan;
        sumN(L)=nan;

        %%%%%%%%
        %OUTPUT:
        %%%%%%%%
        Ex  = Am;
        Ey  = Bm;
        Exy = ABm;
        Vx  = Av;
        Vy  = Bv;
        N   = sumN;
    
    elseif strcmp(lower(mode(1:7)),'mensual') %CORRELACION MENSUAL

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %6. CALCULO LAS ESPERANZAS (MEDIA Y VARIANZA):
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %======================
        %1. OBTENGO LAS MEDIAS:
        %======================
        Am=ones(m,n,p)*nan;
        Bm=ones(m,n,p)*nan;
        A2m=ones(m,n,p)*nan;
        B2m=ones(m,n,p)*nan;
        ABm=ones(m,n,p)*nan;
        
        %ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
        %COJO SOLO LAS POS(i,j,k) DONDE SE TIENEN UN MINIMO DE VALORES
        %PREDEFINIDO, PARA EL CALCULO DE LAS CORRS. MENSUALES:
        %ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
        I=find(N>=0.5*h*h); %33%(o sea, 16 pixels del 7x7=49) con datos.

        Am(I)=CA0(I)./N(I);        %E(A) %Am(55:65,200:210)
        Bm(I)=CB0(I)./N(I);        %E(B)
        A2m(I)=CA0A0(I)./N(I);     %E(A2)
        B2m(I)=CB0B0(I)./N(I);     %E(B2)
        ABm(I)=CA0B0(I)./N(I);     %E(AB)
        
        %=============================================
        %2. OBTENGO LAS VARIANZAS: VAR(X)=E(X2)-E2(X):
        %=============================================
        Av=A2m-(Am.^2);    %Var(A)
        Bv=B2m-(Bm.^2);    %Var(B)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PONGO LA MASCARA TERRESTRE:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L3D=[];
        for k=1:p
            L3D=[L3D;L+(m*n)*(k-1)];
        end
        Am(L3D)=nan;
        Bm(L3D)=nan;
        ABm(L3D)=nan;
        Av(L3D)=nan;
        Bv(L3D)=nan;
        N(L3D)=nan;
        %??????????????????????????????????????????????????????????????????
% $$$         M=zeros(180,360,12);
% $$$         M(L3D)=1;
% $$$         for k=1:12,subplot(3,4,k),imagesc(M(:,:,k)),colorbar('horiz'),end
        %??????????????????????????????????????????????????????????????????

        %%%%%%%%
        %OUTPUT:
        %%%%%%%%
        Ex  = Am;
        Ey  = Bm;
        Exy = ABm;
        Vx  = Av;
        Vy  = Bv;
    
    end %fin corr.anual/mensual
    
else %X e Y son matrices 2D
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %6. CALCULO LAS ESPERANZAS (MEDIA Y VARIANZA):
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %======================
    %1. OBTENGO LAS MEDIAS:
    %======================
    Am=ones(m,n)*nan;
    Bm=ones(m,n)*nan;
    A2m=ones(m,n)*nan;
    B2m=ones(m,n)*nan;
    ABm=ones(m,n)*nan;
    
    %ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
    %COJO SOLO LAS POS(i,j,k) DONDE SE TIENEN UN MINIMO DE VALORES
    %PREDEFINIDO, PARA EL CALCULO DE LAS CORRS. MENSUALES:
    %ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
%    I=find(N>=0.5*h*h); %33%(o sea, 16 pixels del 7x7=49) con datos.
    I=find(N>=1);
    
    Am(I)=CA0(I)./N(I);        %E(A) %Am(55:65,200:210)
    Bm(I)=CB0(I)./N(I);        %E(B)
    A2m(I)=CA0A0(I)./N(I);     %E(A2)
    B2m(I)=CB0B0(I)./N(I);     %E(B2)
    ABm(I)=CA0B0(I)./N(I);     %E(AB)
    
    %=============================================
    %2. OBTENGO LAS VARIANZAS: VAR(X)=E(X2)-E2(X):
    %=============================================
    Av=A2m-(Am.^2);    %Var(A)
    Bv=B2m-(Bm.^2);    %Var(B)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PONGO LA MASCARA TERRESTRE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Am(L)=nan;
    Bm(L)=nan;
    ABm(L)=nan;
    Av(L)=nan;
    Bv(L)=nan;
    N(L)=nan;
    %??????????????????????????????????????????????????????????????????
% $$$     length(L)
% $$$     M=zeros(180,360);
% $$$     M(L)=1;
% $$$     imagesc(M)
% $$$     colorbar('horiz')
% $$$     pause
    %??????????????????????????????????????????????????????????????????

    %%%%%%%%
    %OUTPUT:
    %%%%%%%%
    Ex  = Am;
    Ey  = Bm;
    Exy = ABm;
    Vx  = Av;
    Vy  = Bv;

end %fin matrices 3D/2D

%%%%%%%%%%
%GRAFICAS:
%%%%%%%%%%
% $$$ MAP=jet;
% $$$ MAP(1,:)=[0 0 0];
% $$$ figure(1)
% $$$ imagesc(Ex)
% $$$ colormap(MAP)
% $$$ colorbar('horiz')
% $$$ 
% $$$ figure(2)
% $$$ imagesc(Ey)
% $$$ colormap(MAP)
% $$$ colorbar('horiz')
% $$$ 
% $$$ figure(3)
% $$$ imagesc(Vx)
% $$$ colormap(MAP)
% $$$ colorbar('horiz')
% $$$ 
% $$$ figure(4)
% $$$ imagesc(Vy)
% $$$ colormap(MAP)
% $$$ colorbar('horiz')

display('fin corr_estadisticos')
