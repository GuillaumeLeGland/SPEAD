function [keyKN,keyTOPT,key2T,keyModelResol,keyPhysics,keySinking,keyPARseasonality,keyPARextinction,keyNutrientSupplyFlux,keyNutrientPulses,keyKTW,keyTraitDiffusion,keyGalfaConstant,keyTraitDiffusionConstant,keyNutrientSupplyFrequencyConstant,keyAssimConstant] = jamstecrest_gaussecomodel1D_keys()

%========================================================================
% Keys to activate different physical processes or numerical schemes
%........................................................................
%%keyTraitAxis = 'ESD';
%%keyTraitAxis = 'SST';
% Choose 1 trait, 2 traits, or both (Le Gland, 22/11/2019)
keyKN   = 'yes'; % 1-trait model with KN as a trait. No temperature effect.
keyTOPT = 'yes'; % 1-trait model with TOPT as a trait. Standard KN is chosen 
key2T   = 'yes'; % 2-trait model
%........................................................................
%%keyModelResol='0D';
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
keyNutrientSupplyFlux = 'not'; % If you want nutrient supply (and dilution)
keyNutrientPulses = 'Continuous'; % Added to the file of keys (Le Gland, 28/06/2019)
%keyNutrientPulses = 'Discrete';
%........................................................................
keyKTW = 'not';
keyTraitDiffusion = 'yes';
%........................................................................
keyTraitDiffusionConstant = 'yes';
%........................................................................
keyGalfaConstant = 'yes';
%........................................................................
keyNutrientSupplyFrequencyConstant = 'yes'; %If you do *not* want depth-increasing pulse frequency. 
%........................................................................
keyAssimConstant = 'yes';
%........................................................................
%========================================================================
%************************************************************************
return


