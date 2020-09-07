function [d1uxdxSym,d2uxdxSym,d3uxdxSym,d4uxdxSym] = jamstecrest_gaussecomodel1D_symbolicderivatives()
%========================================================================
%........................................................................
%DEFINE SYMBILIC VARIABLES:
syms mup  knp  amup aknp bmup 
syms muxj knxj qxj uxj xj 
syms din 
%........................................................................
%%muxj = mup .* exp(amup*xj); %Phy maximum grazing rate as a function of cell size [d-1] 
muxj = mup .* exp(amup*xj + bmup*xj^2);
knxj = knp .* exp(aknp*xj); %Phy half-sat uptake as a function of cell size [mmolN*m-3] 
%........................................................................
d1muxdxSym = diff(muxj,xj,1); %1ST derivative.
d2muxdxSym = diff(muxj,xj,2); %2ND derivative. 
d3muxdxSym = diff(muxj,xj,3); %3RD derivative.
d4muxdxSym = diff(muxj,xj,4); %4TH derivative.
%........................................................................
qxj = (din ./ (din + knxj)); 
uxj = muxj .* qxj; 
%........................................................................
d1uxdxSym = diff(uxj,xj,1); %1ST derivative.
d2uxdxSym = diff(uxj,xj,2); %2ND derivative. 
d3uxdxSym = diff(uxj,xj,3); %3RD derivative.
d4uxdxSym = diff(uxj,xj,4); %4TH derivative.
%........................................................................
%========================================================================
%************************************************************************
return

