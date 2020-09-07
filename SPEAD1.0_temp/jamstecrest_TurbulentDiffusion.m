%function [DIFF] = jamstecrest_TurbulentDiffusion(Cz,dt,dz,kz,kzI,nz,keyBoundaryCondition)
%kzI is no longer necessary (Le Gland, 10/09/2019)
function [DIFF] = jamstecrest_TurbulentDiffusion(Cz,dt,dz,kz,nz,keyBoundaryCondition)
%function [DIFF] = jamstecrest_TurbulentDiffusion(Cz,dt,dz,kz,nz,split,keyScheme)

%Add implicit mode to avoid instability, as in NEMO 
%When time step is much larger than dz^2/(2*KZ), the explicit scheme is
%unstable, the Crank-Nicolson scheme creates oscillations and the implicit
%scheme homgenizes the concentrations. This is why we need implicit.
%Reducing the time step if it only affects vertical diffusion but not the
%(slower) ecological dynamics is not relevant. (Le Gland, 15/10/2019)
% Link between implicit annr Runge Kutta ?
% For now, only works with Neumann boundary condition

% if strcmp(keyScheme,'Implicit')

T = zeros(nz); % Inverse of the transport matrix (Le Gland, 15/10/2019)
for i=2:nz-1
    T(i-1,i) = - kz(i) * dt / dz^2;
    T(i+1,i) = - kz(i+1) * dt / dz^2;
    T(i,i) = 1 - T(i-1,i) - T(i+1,i);
end
T(2,1) = - kz(2) * dt / dz^2;
T(1,1) = 1 - T(2,1);
T(nz-1,nz) = - kz(nz) * dt / dz^2;
T(nz,nz) = 1 - T(nz-1,nz);

DIFF = (T \ Cz - Cz) / dt;

% end

% Explicit reducing the time step only for vertical diffusion is another
% option (Le Gland, 16/10/2019). It is not even more time-expansive !
% The difference with the implicit scheme is insignificant !
% if strcmp(keyScheme,'Explicit')
% split = 500;
% for i=2:nz-1
%     T(i-1,i) = kz(i) * dt / (split*dz^2);
%     T(i+1,i) = kz(i+1) * dt / (split*dz^2);
%     T(i,i) = 1 - T(i-1,i) - T(i+1,i);
% end
% T(2,1) = kz(2) * dt / (split*dz^2);
% T(1,1) = 1 - T(2,1);
% T(nz-1,nz) = kz(nz) * dt / (split*dz^2);
% T(nz,nz) = 1 - T(nz-1,nz);
% 
% DIFF = ((T^split) * Cz - Cz) / dt;
% end

return

Cz  = Cz(:)'; %row vector (0-Zm).
kz  = kz(:)'; %row vector (0-Zm).
%kzI = kzI(:)'; %row vector (0-Zm).



%%%%%%%%%%%%%%%%%%%%%%%
%FOR THE SURFACE (j=1):
%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
J = [1];
%...................................................................................
if strcmp(keyBoundaryCondition,'Reflectante')
    Czsurface = Cz(J); %Frontier Boundary Condition (reflectance frontier - Neumann)
elseif strcmp(keyBoundaryCondition,'Absorbente')
    Czsurface = 0; %Frontier Boundary Condition (absorbent frontier - Dirichlet)
end
%...................................................................................
kzsurface = kz(J);
%...................................................................................
%DIFFsurface = (kzI(J)*(Cz(J+1)-Cz(J)) - kzsurface*(Cz(J)-Czsurface)) / dz^2; %kzI(1) = kz(j=1/2)
DIFFsurface = (kz(J+1)*(Cz(J+1)-Cz(J)) - kzsurface*(Cz(J)-Czsurface)) / dz^2;
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MIDDLE DEPTHS [j = 2:nz-1]:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
J = [2:nz-1];
%...................................................................................
%DIFFmiddle = (kzI(J).*(Cz(J+1)-Cz(J)) - kzI(J-1).*(Cz(J)-Cz(J-1))) / dz^2;
DIFFmiddle = (kz(J+1).*(Cz(J+1)-Cz(J)) - kz(J).*(Cz(J)-Cz(J-1))) / dz^2;
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%
%BOTTOM (j=nz):
%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
J = [nz];
%...................................................................................
if strcmp(keyBoundaryCondition,'Reflectante')
    Czdeep = Cz(J); %Frontier Boundary Condition (reflectance frontier - Neumann)
elseif strcmp(keyBoundaryCondition,'Absorbente')
    Czdeep = 0; %Frontier Boundary Condition (absorbent frontier - Dirichlet)
end
%...................................................................................
%kzdeep = kz(J);
kzdeep = kz(J+1);
%...................................................................................
%DIFFdeep = (kzdeep*(Czdeep-Cz(J)) - kzI(J-1)*(Cz(J)-Cz(J-1))) / dz^2;
DIFFdeep = (kzdeep*(Czdeep-Cz(J)) - kz(J)*(Cz(J)-Cz(J-1))) / dz^2;
%...................................................................................
%===================================================================================

%%%%%%%%%%%%%%%
%TOTAL PROFILE:
%%%%%%%%%%%%%%%
DIFF = [DIFFsurface,DIFFmiddle,DIFFdeep]; %[mmolN*m-3*d-1]

%%%%%%%%
%OUTPUT:
%%%%%%%%
DIFF = DIFF(:); %column vector.

