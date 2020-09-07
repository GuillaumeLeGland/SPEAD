function [jac]=mynumjacobian(f,x)
%------------------------------------------------
% Compute a forward difference Jacobian f'(x).
%
% [jac]=mynumjacobian(f,x)
%
% This code comes with no guarantee or warranty of any kind.
%
% inputs:
%         x  = steady state point.
%         f  = ODE function (name in strings)
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

n=length(x);
%...............
% $$$ epsnew=sqrt(eps); %difference increment.
% $$$ epsnew=1d-7; 
epsnew=1d-4; %USAR ESTE.
%...............

f0=feval(f,0,x); %where "f" is called with ode45.
for j=1:n
    w=zeros(n,1);
    w(j)=1; %direction.

    %Scale the step:
    epsnew = epsnew/norm(w);
    if norm(x) > 0
	epsnew=epsnew*norm(x);
    end
    del=x+epsnew*w;
    f1=feval(f,0,del); %where "f" is called with ode45.
    z = (f1 - f0)/epsnew;
    jac(:,j)=z;
end
