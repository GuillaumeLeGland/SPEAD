function [d1uxdxSym,d2uxdxSym,d3uxdxSym,d4uxdxSym] = jamstecrest_derivatives_symbolic()

%DEFINE SYMBILIC VARIABLES:
syms mu0 ks0 amu aks 
syms mux ksx qx ux x 
syms DIN 

mux = mu0 * exp(amu*x); 
ksx = ks0 * exp(aks*x); 

qx = (DIN ./ (DIN + ksx)); 
ux = mux .* qx; 

d1uxdxSym = diff(ux,x,1); 
d2uxdxSym = diff(ux,x,2); 
d3uxdxSym = diff(ux,x,3); 
d4uxdxSym = diff(ux,x,4); 

return

