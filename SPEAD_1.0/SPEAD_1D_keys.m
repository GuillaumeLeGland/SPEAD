function [keyKN,keyTOPT,key2T,keyDisc,keyModelResol,keyPhysics,keySinking,keyPARseasonality,keyPARextinction,keyKTW,keyTraitDiffusion] = SPEAD_1D_keys()

%========================================================================
% Keys to activate different physical processes or numerical schemes
%........................................................................
% Run the single-trait models, the double trait model or all of them
keyKN   = 'yes'; % 1-trait aggregate model with KN as a trait. Optimal temperature is equal to environment temperature.
keyTOPT = 'yes'; % 1-trait aggregate model with TOPT as a trait. Kn is equal to the environment DIN. 
key2T   = 'yes'; % 2-trait aggregate model
keyDisc = 'yes'; % 2-trait discrete model
%........................................................................
%keyModelResol = '0D'; % The model runs but the analysis and figures will not work
keyModelResol='1D'; %Depth-resolved (several nodes in depth).
%........................................................................
if strcmp(keyModelResol,'0D')
    keyPhysics='not'; % Never change this one
else
    keyPhysics='yes'; %If you want physical processes (turbulent diffusion).
end
%........................................................................
keySinking='yes'; %If you want advective processes (vertical sinking).
%........................................................................
keyPARseasonality='yes'; %If you want seasonal solar radiation (PAR).
%........................................................................
keyPARextinction='yes'; %If you want depth-decreasing solar radiation (PAR).
%........................................................................
keyKTW = 'not';
keyTraitDiffusion = 'yes';
%........................................................................
%========================================================================
%************************************************************************
return


