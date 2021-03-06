function [varargout]=myEasyFit(xinput,yinput)

%==========================================
% $$$ %..................
% $$$ xmax = max(xinput);
% $$$ ymax = max(yinput);
% $$$ %..................
% $$$ x = xinput;
% $$$ y = yinput;
% $$$ %..................
% $$$ plot(x,y)
% $$$ %..................
% $$$ H = showfit('a*(x/(b + x))')
% $$$ %..................
% $$$ a = H.m(1);
% $$$ b = H.m(2); 
% $$$ %..................
%==========================================
% $$$ %..................
% $$$ xmax = max(xinput);
% $$$ ymax = max(yinput);
% $$$ %..................
% $$$ x = xinput./xmax;
% $$$ y = yinput./ymax;
% $$$ %..................
% $$$ plot(x,y)
% $$$ %..................
% $$$ H = showfit('x/(b + x)')
% $$$ %..................
% $$$ bb = H.m; 
% $$$ %..................
% $$$ a = ymax;
% $$$ b = bb*xmax
% $$$ %..................
%==========================================
%..................
xmax = max(xinput);
ymax = max(yinput);
%..................
x = xinput./xmax;
y = yinput./ymax;
%..................
plot(x,y)
%..................
H = showfit('a*(x/(b + x))')
%..................
aa = H.m(1); 
bb = H.m(2); 
%..................
a = aa*ymax;
b = bb*xmax
%..................
rho = H.r;
r2  = rho^2;
%..................
%==========================================
varargout{1}=a;
varargout{2}=b;
varargout{3}=r2;
%==========================================

