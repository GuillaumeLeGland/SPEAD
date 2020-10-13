function [Ndiff] = SPEAD_discretemodel1D_traitDispersion(Nxy,rxy,nux,nuy,nxphy,nyphy,dx,dy) 
%============================================================================
%DEFINE TRAIT AXIS MIN-MAX VALUES:
%............................................................................
imin = 1; 
imax = nxphy;
jmin = 1;
jmax = nyphy;
%............................................................................
%============================================================================
N_ex  = [Nxy(imin,:);Nxy;Nxy(imax,:)];
N_exy = [N_ex(:,jmin),N_ex,N_ex(:,jmax)]; % Population concentration [mmolN * m-3]
r_ex  = [rxy(imin,:);rxy;rxy(imax,:)];
r_exy = [r_ex(:,jmin),r_ex,r_ex(:,jmax)]; % Growth rate [d-1]
%............................................................................
Nr_xy = N_exy .* r_exy; % Population growth [mmolN * m-3 * d-1]
%............................................................................
%============================================================================
%----------------------------------------------------------------------------
%Delta production units: ([d-1] x [mmolN * m-3]) / [trait-1] 
%----------------------------------------------------------------------------
delta_production_x1 = (Nr_xy(1:end-2,2:end-1) - Nr_xy(2:end-1,2:end-1)) / dx; %[mmolN * m-3 * d-1] x [trait-1] 
delta_production_x2 = (Nr_xy(2:end-1,2:end-1) - Nr_xy(3:end,2:end-1)) / dx; %[mmolN * m-3 * d-1] x [trait-1]
delta_production_y1 = (Nr_xy(2:end-1,1:end-2) - Nr_xy(2:end-1,2:end-1)) / dy; %[mmolN * m-3 * d-1] x [trait-1] 
delta_production_y2 = (Nr_xy(2:end-1,2:end-1) - Nr_xy(2:end-1,3:end)) / dy; %[mmolN * m-3 * d-1] x [trait-1]
%............................................................................
%----------------------------------------------------------------------------
% Flux units: [trait^2] x ([mmolN * m-3 * d-1] x [trait-1]) x [trait-1] = [mmolN * m-3 * d-1]
%----------------------------------------------------------------------------
% Note that here I could impose a different nu for x and y
% This is done since 10/09/2019
Flux_x1 = (nux * delta_production_x1) / dx; %Source [mmolN * m-3 * d-1]
Flux_x2 = (nux * delta_production_x2) / dx; %Sinks [mmolN * m-3 * d-1]
Flux_y1 = (nuy * delta_production_y1) / dy; %Source [mmolN * m-3 * d-1]
Flux_y2 = (nuy * delta_production_y2) / dy; %Sinks [mmolN * m-3 * d-1]
%............................................................................
%============================================================================
%COMPUTE SOURCE-SINKS BALANCE FOR EACH GRID POINT:
%............................................................................
Ndiff = (Flux_x1 - Flux_x2) + (Flux_y1 - Flux_y2); %[mmolN * m-3 * d-1]

%............................................................................
%============================================================================
%****************************************************************************
return

