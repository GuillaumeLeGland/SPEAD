function [RP,RPsig,TP,LtotP]=corr_rpearson(Ex,Ey,Exy,Vx,Vy,N,L,h);
%-------------------------------------------------------
%NOTA: (pag. 583. SOKAL)
%1. R-PEARSON: % r=(E(xy)-E(x)E(y))/sqrt(Var(x)*Var(y))    
%2. SIGNIFICANCIA: (Test "ts": Ho: nu=0 ; H1: nu~=0;)
%Sr = sqrt((1-r^2)/(n-2)); 
%ts = (r - nu)/Sr = (r - 0)/Sr = r*sqrt((n-2)/(1-r^2));
%-------------------------------------------------------    
[m,n,p]=size(Ex);

if p==1 %Ex es una matriz 2D
    R = ones(180,360)*nan;
    T = ones(180,360)*nan;
    if h==3
        J=find(N < 0.5*h*h*p);
    elseif h>3
        J=find(N < 100); %donde haya menos de 100 ptos no calculo correlacion.
    end
elseif p>1 %Ex es una matriz 3D.
    R = ones(180,360,p)*nan;
    T = ones(180,360,p)*nan;
    J=find(N < 0.5*h*h);
end
    
R = (Exy - Ex.*Ey) ./ sqrt(Vx.*Vy);
T = R.*sqrt((N-2)./(1-R.^2));

%======================================
%PONGO NAN DONDE R ES NO SIGNIFICATIVA:
%======================================
Rsig=R;
Inosig=find(T>-1.9 & T<+1.9);
Rsig(Inosig)=nan;

%===============================================
%ELIMINO LAS POSICIONES QUE SE REPITEN EN L y J:
%===============================================
if p==1
    M=zeros(180,360);
    M([L;J])=1; %Land + No-enough-data = 1.
    Ltot=find(M==1); %vetor-col. con las pos. en la matriz 2D "RP" con
                     %Land y no-enough-data.
elseif p>1
    Ltot=ones(180*360,p)*nan;
    M=zeros(180,360,p);
    M([L;J])=1;
    for k=1:12
	F=find(M(:,:,k)==1);
	Ltot(F,k)=F; %12 vectores-col con las pos. en la matriz 3D "RP"
                     %con Land y no-enough-data (cada vector corresponde
                     %a cada mes).
    end
end

%OUTPUT:
RP=R;
RPsig=Rsig;
TP=T;
LtotP=Ltot;

%***************************
display('fin corr_rpearson.m')
