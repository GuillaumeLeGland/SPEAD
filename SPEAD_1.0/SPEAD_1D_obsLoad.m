function [CHL_obs,PP_obs,NO3_obs,PON_obs] = SPEAD_1D_obsLoad( )
%--------------------------------------------------------------------------
% Load all BATS observations to be compared with model outputs

PP_obs  = load('INPUTS/PPbatsClimMonthlyProfilesPolyregress.txt');
CHL_obs = load('INPUTS/CHLbatsClimMonthlyProfilesPolyregress.txt');
NO3_obs = load('INPUTS/NO3batsClimMonthlyProfilesPolyregress.txt');
PON_obs = load('INPUTS/PONbatsClimMonthlyProfilesPolyregress.txt');
PP_obs  = PP_obs * (1/12) * (16/106); % Converts mgC/m3/day to mmol/m2/day
PON_obs = PON_obs / 14; % Converts ug/kg to mmol/m3

end

