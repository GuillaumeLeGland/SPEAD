%function [Ndiff] = jamstecrest_discretemodel1D_traitDispersion(Nx,rx,nux,nphy,dx)
% Case with 2 traits (Le Gland, 16/07/2019)
%function [Ndiff] = jamstecrest_discretemodel1D_traitDispersion(Nxy,rxy,nuxy,nxphy,nyphy,dx,dy) 
function [Ndiff] = jamstecrest_discretemodel1D_traitDispersion(Nxy,rxy,nux,nuy,nxphy,nyphy,dx,dy) 
%%function [Ndiff] = evobranchAD_randomizedDispers1D(Nx,rx,nux,nphy,dx)

%============================================================================
%DEFINE TRAIT AXIS MIN-MAX VALUES:
%............................................................................
imin = 1; 
%imax = nphy; 
% Case with 2 traits (Le Gland, 16/07/2019)
imax = nxphy;
jmin = 1;
jmax = nyphy;
%............................................................................
%============================================================================
%COMPUTE BIOMASS FLUXES AMONG GRID POINTS:
%............................................................................
%Extend domain with two grid-points by repeating boundary conditions: 
%N_ex = [Nx(imin);Nx(:);Nx(imax)]; %Population concentration [mmolN * m-3] 
%r_ex = [rx(imin);rx(:);rx(imax)]; %Growth rate [d-1] 
%............................................................................
%N_1 = N_ex(1:end-2); %Biomass at N(i-1)
%N_2 = N_ex(2:end-1); %Biomass at N(i) -- middle point
%N_3 = N_ex(3:end);   %Biomass at N(i+1)
%............................................................................
%r_1 = r_ex(1:end-2); %Specific growth rate at rx(i-1)
%r_2 = r_ex(2:end-1); %Specific growth rate at rx(i) -- middle point
%r_3 = r_ex(3:end);   %Specific growth rate at rx(i+1)
%............................................................................
%============================================================================
%----------------------------------------------------------------------------
%Delta production units: ([d-1] x [mmolN * m-3]) / [trait-1] 
%----------------------------------------------------------------------------
%delta_production_1 = ((r_1.*N_1) - (r_2.*N_2)) / dx; %[mmolN * m-3 * d-1] x [trait-1] 
%delta_production_2 = ((r_2.*N_2) - (r_3.*N_3)) / dx; %[mmolN * m-3 * d-1] x [trait-1] 
%............................................................................
%----------------------------------------------------------------------------
% Flux units: [trait^2] x ([mmolN * m-3 * d-1] x [trait-1]) x [trait-1] = [mmolN * m-3 * d-1]
%----------------------------------------------------------------------------
%Flux_1 = (nux * delta_production_1) / dx; %Source [mmolN * m-3 * d-1]
%Flux_2 = (nux * delta_production_2) / dx; %Sinkes [mmolN * m-3 * d-1]
%............................................................................
%============================================================================
%COMPUTE SOURCE-SINKS BALANCE FOR EACH GRID POINT:
%............................................................................
%Ndiff = (Flux_1 - Flux_2); %[mmolN * m-3 * d-1]
%............................................................................
%% $$$ %%Ndiff_bis = nux * ((r_1.*N_1) - (r_3.*N_3)); %[mmolN * m-3 * d-1] %WRONG!!!!
%% $$$ %............................................................................
%% $$$ Ndiff_bis = nux * ((r_1.*N_1) - 2*(r_2.*N_2) + (r_3.*N_3)); %[mmolN * m-3 * d-1] %Okay 
%% $$$ %............................................................................
%% $$$ show_Ndiff = [Ndiff,Ndiff_bis]
%% $$$ pause

% Case with 2 traits (Le Gland, 16/07/2019)
% Extend domain on each side
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

