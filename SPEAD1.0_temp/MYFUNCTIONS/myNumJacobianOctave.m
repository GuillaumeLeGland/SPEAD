function [jac]=myNumJacobianOctave(f,x)
global keyOctaveMatlabASK
%------------------------------------------------
% Compute a forward difference Jacobian f'(x).
%
% [jac]=myNumJacobian(f,x)
%
% This code comes with no guarantee or warranty of any kind.
%
% inputs:
%         x  = steady state point.
%         f  = ODE function (name in strings) which is called with ode45.
%
% output: 
%         jac = numerical jacobian.
%
% Uses finite difference directional derivative
% Approximate f'(x) w
% 
% Sergio M. Vallina, Dec 10, 2008.
% (modified from C. T. Kelley, November 25, 1993)
%------------------------------------------------

x=x(:);
n=length(x);
%...............
% $$$ epsnew=sqrt(eps); %difference increment. %Original.
epsnew=1d-7; 
% $$$ epsnew=1d-4; %USAR ESTE (if necessary).
%...............
if strcmp(keyOctaveMatlabASK,'Octave')
    f0=feval(f,x,0); %where "f" is called with lsode.
elseif strcmp(keyOctaveMatlabASK,'Matlab')
    f0=feval(f,0,x); %where "f" is called with ode45.
end
%...............
for j=1:n
    %........................
    w=zeros(n,1);
    w(j)=1; %direction.
    %........................
    %Scale the step:
    epsnew = epsnew/norm(w);
    if norm(x) > 0
	epsnew=epsnew*norm(x);
    end
    %........................
    xdel = x + epsnew*w;
    %........................
    if strcmp(keyOctaveMatlabASK,'Octave')
	f1=feval(f,xdel,0); %where "f" is called with lsode.
    elseif strcmp(keyOctaveMatlabASK,'Matlab')
	f1=feval(f,0,xdel); %where "f" is called with ode45.
    end
    %........................
    z = (f1 - f0)/epsnew;
    jac(:,j)=z;
    %........................
end
