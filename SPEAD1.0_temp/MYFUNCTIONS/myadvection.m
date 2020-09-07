function [ADV]=myadvection(C0,dz,w,zjmax)
%...
C0=C0(:)'; %vector fila (0-Zm).
%...
    
CC0=[0,C0,0]; %Anado dos C.F.=0
J=[1:zjmax];
JJ=J+1;%[2:zjmax+1]; %como JJ se aplicara a CC0 (que tiene 2 ptos max como CF) necesito desplazar los subindices un posicion para que CCO(JJ)=C0(J).
%..............................................
conc0=C0(J);
conc0bis=CC0(JJ);
if conc0~=conc0bis
    error('C0(J) y CC0(JJ) deben ser iguales!')
end
%..............................................
if w>0 %sinking.
    %ALL NODES EXCEPT FRONTERA:
    entra=w*CC0(JJ-1);
    sale=w*CC0(JJ);
    %FRONTERA (deep): 
    %a) Si quiero que sea reflectante pongo sale(end)=0.
    %b) Si quiero que sea absorbente pongo sale(end)=sale(end) (es decir no cambio nada...).
    sale(end)=0;
%    sale(end)=sale(end);
elseif w<0 %upwelling.
    entra=-w*CC0(JJ+1);
    sale=-w*CC0(JJ);
    %FRONTERA (surface): 
    %a) Si quiero que sea reflectante pongo sale(1)=0.
    %b) Si quiero que sea absorbente pongo sale(1)=sale(1) (es decir no cambio nada...).
%    sale(1)=0;
    sale(1)=sale(1);
end
ADV=(entra-sale)./dz;
ADV=ADV(:); %vector columna.
