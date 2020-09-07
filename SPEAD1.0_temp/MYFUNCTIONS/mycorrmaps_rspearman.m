function [RS,RSsig,TS]=mycorrmaps_rspearman(VARXm,VARYm)
   
[m,n,p]=size(VARXm);
if p==1
    error('el VARXm debe ser una matriz cubica')
end

display('start mycorrmaps_rspearman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%7. OBTENGO "EL RANKING" (DE MENOR A MAYOR) DE LOS VALORES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------
%Vi  : Valores del vector ordenados de menor a mayor.
%RGi : Posiciones originales que ocupaban en el vector antes de ordenarse.
%-------------------------------------------------------------------------
[V1,E1]=sort(VARXm,3); %E1 = Etiquetas (pos.orig. de los valores en los vect de VARXm)
[V2,E2]=sort(VARYm,3);
[X1,RG1]=sort(E1,3); %obtengo los rangos (ranks) (la X1 no se usa)
[X2,RG2]=sort(E2,3);
clear X1 X2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO LOS GRADOS DE LIBERTAD (No DE PUNTOS):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=zeros(180,360,p);
B=zeros(180,360,p);

IA=find(VARXm>=0);
IB=find(VARYm>=0);

if IA==IB
    I=IA;
else
   error('deben ser iguales')
end
A(I)=1; %pongo 1's donde hay valores en Am (el resto es 0's)
B(I)=1;

NA = sum(A,3);
NB = sum(B,3);
if NA==NB
    N=NA;
else
    error('Error!: tienen que ser iguales')
end

%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
%BUSCO LAS POS(lat,long) "J" DONDE HAYA MENOS DE 6 MESES EN LOS
%QUE SE HA PODIDO CALCULAR LA MEDIA ESPACIAL. EN ELLAS APLICARE UNA
%MASCARA NEGRA A LA "R Spearman".
%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
%jmin=6;
jmin=7; %prefiero que haya al menos 7 meses con dato.
J=find(N<jmin);
N(J)=nan; %solo me quedo con los vectores con mas de 6 medias.

%æææææææææææææææææææ
%CALCULO R-Spearman:
%æææææææææææææææææææ
D=(RG1-RG2).^2;
R = 1 - ( 6*sum(D,3) ./ (N.*(N.^2-1)) ); %correlacion spearman 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%8. SIGNIFICANCIA: (Test "ts": Ho: nu=0 ; H1: nu~=0;)
%    Sr = sqrt((1-r^2)/(n-2)); 
%    ts = (r - nu)/Sr = (r - 0)/Sr = r*sqrt((n-2)/(1-r^2));
%NOTA: Se aplica el mismo test que para el coef.corr.PEARSON (parametrico)
%aunque estemos con el de SPEARMAN (no-param.) (ver: pag.607. SOKAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = R.*sqrt((N-2)./(1-R.^2));

%======================================
%PONGO NAN DONDE R ES NO SIGNIFICATIVA:
%======================================
Inosig=find(T>-1.9 & T<+1.9); %t=1.9 implica alfa=0.05 (see SOKAL, pag:102).
Rsig=R;
Rsig(Inosig)=nan; 

%%%%%%%%
%OUTPUT:
%%%%%%%%
RS=R;
RSsig=Rsig;
TS=T;

%æææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
%PONGO EN PIXEL SIN DATO EL VALOR -1 Y OTRO CON VALOR +1 PARA OBTENER
%LUEGO UNA COLORBAR GUAY:
%æææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
RS(1,1)=-1;
RS(1,2)=+1;
RSsig(1,1)=-1;
RSsig(1,2)=+1;
