function [ticks] = mycolorbarticks(Adata,nTicks,nDps)

% Find the limits
lims = [min(Adata(:)) max(Adata(:))]; 

% Function for rounding to specified decimal places
dprnd = @(x, dps)round(x*(10.^dps))./(10.^dps);

% Generate ticks
% $$$ nDps = 2; %Decimal precision. 
% $$$ nTicks = 5; %Number of ticks.
ticks = dprnd(linspace(lims(1), lims(2), nTicks), nDps);
