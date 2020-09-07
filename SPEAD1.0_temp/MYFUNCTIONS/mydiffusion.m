function [DIFF]=mydiffusion(C0,dz,kz,kzI,zjmax)

C0=C0(:)'; %vector fila (0-Zm).
kz=kz(:)'; %vector fila (0-Zm).
kzI=kzI(:)'; %vector fila (0-Zm).
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULO TURBULENT DIFFUSION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%æææææææææææææææææææææææææ
%PARA LA SUPERFICIE (j=1):
%æææææææææææææææææææææææææ
J=1;
C0surface=C0(J); %C.F. (reflectante)
%C0surface=0;
kzsurface=kz(J);
DIFFsurface=(kzI(J)*(C0(J+1)-C0(J))-kzsurface*(C0(J)-C0surface))/dz^2; %kzI(1) = kz(j=1/2)

%ææææææææææææææææææææææææææææææææææææææ
%PARA LAS DEPTHS EN EL MEDIO (j=2:zjmax-1):
%ææææææææææææææææææææææææææææææææææææææ
J=2:zjmax-1;
DIFFmedio=(kzI(J).*(C0(J+1)-C0(J))-kzI(J-1).*(C0(J)-C0(J-1)))/dz^2;

%ææææææææææææææææææææææ
%PARA EL FONDO (j=zjmax):
%ææææææææææææææææææææææ
J=zjmax;
C0deep=C0(J); %C.F.
%C0deep=0; %C.F.
kzdeep=kz(J);
DIFFdeep=(kzdeep*(C0deep-C0(J))-kzI(J-1)*(C0(J)-C0(J-1)))/dz^2;

%æææææææææææææ
%PERFIL TOTAL:
%æææææææææææææ
DIFF=[DIFFsurface,DIFFmedio,DIFFdeep];
DIFF=DIFF(:); %vector columna.
