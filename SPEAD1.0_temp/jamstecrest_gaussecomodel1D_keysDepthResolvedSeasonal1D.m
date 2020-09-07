function [keyTraitAxis,keyModelResol,keyFastNumericalSolving,keyPhysics,keySinking,keyPARseasonality,keyPARextinction,keyNutrientSupplyFlux,keyNutrientPulses,keyKTW,keyTraitDiffusion,keyGalfaConstant,keyTraitDiffusionConstant,keyLogBase,keyNutrientSupplyFrequencyConstant] = jamstecrest_gaussecomodel1D_keysDepthResolvedSeasonal1D()

%========================================================================
% Keys to activate different physical processes or numerical schemes
%........................................................................
%%keyTraitAxis = 'ESD';
keyTraitAxis = 'SST';
%........................................................................
%%keyLogBase = 'Base2'; 
keyLogBase = 'BaseExp'; %(I THINK IT NEEDS TO BE LOGEXP FOR THE ANALYTIC APPROACH TO BE VALID).
%........................................................................
%keyModelResol='0D';
keyModelResol='1D'; %Depth-resolved (several nodes in depth).
%........................................................................
keyFastNumericalSolving='yes';
%........................................................................
keyPhysics='yes'; %If you want physical processes (turbulent diffusion).
%........................................................................
keySinking='yes'; %If you want advective processes (vertical sinking).
%........................................................................
keyPARseasonality='yes'; %If you want seasonal solar radiation (PAR).
%........................................................................
keyPARextinction='yes'; %If you want depth-decreasing solar radiation (PAR).
%........................................................................
keyNutrientSupplyFlux = 'not';
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
%========================================================================
%************************************************************************
return

