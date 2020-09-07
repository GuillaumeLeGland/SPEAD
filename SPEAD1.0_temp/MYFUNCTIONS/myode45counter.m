function [tout, yout] = ode45(F, tspan, y0, tol, trace, count)
% ode45 (v1.0) integrates a system of ordinary differential equations using
% 4th & 5th order embedded Runge-Kutta-Fehlberg formulas.
% This requires 6 function evaluations per integration step.
% This particular implementation uses the 5th order estimate for yout.
%
% The Fehlberg 4(5) coefficients are from a tableu in
% U.M. Ascher, L.R. Petzold, Computer Methods for  Ordinary Differential Equations
% and Differential-Agebraic Equations, Society for Industrial and Applied Mathematics
% (SIAM), Philadelphia, 1998
%
% The error estimate formula and slopes are from
% Numerical Methods for Engineers, 2nd Ed., Chappra & Cannle, McGraw-Hill, 1985
%
% Usage:
%         [tout, yout] = ode45(F, tspan, y0, tol, trace, count)
%
% INPUT:
% F     - String containing name of user-supplied problem description.
%         Call: yprime = fun(t,y) where F = 'fun'.
%         t      - Time (scalar).
%         y      - Solution column-vector.
%         yprime - Returned derivative COLUMN-vector; yprime(i) = dy(i)/dt.
% tspan - [ tstart, tfinal ]
% y0    - Initial value COLUMN-vector.
% tol   - The desired accuracy. (Default: tol = 1.e-6).
% trace - If nonzero, each step is printed. (Default: trace = 0).
% count - if nonzero, variable 'counter' is initalized, made global
%         and counts the number of state-dot function evaluations
%         'counter' is incremented in here, not in the state-dot file
%         simply make 'counter' global in the file that calls ode45
%
% OUTPUT:
% tout  - Returned integration time points (column-vector).
% yout  - Returned solution, one solution column-vector per tout-value.
%
% The result can be displayed by: plot(tout, yout).
%
% Marc Compere
% compere at mail dot utexas dot edu
% 6 October 1999

if nargin < 6, count = 0; end
if nargin < 5, trace = 0; end
if nargin < 4, tol = 1.e-6; end

pow = 1/8;

% The Fehlberg 4(5) coefficients:
 a_(1,1)=0;
 a_(2,1)=1/4;
 a_(3,1)=3/32; 
 a_(3,2)=9/32;
 a_(4,1)=1932/2197; 
 a_(4,2)=-7200/2197; 
 a_(4,3)=7296/2197;
 a_(5,1)=439/216; 
 a_(5,2)=-8; 
 a_(5,3)=3680/513; 
 a_(5,4)=-845/4104; 
 a_(6,1)=-8/27; 
 a_(6,2)=2; 
 a_(6,3)=-3544/2565; 
 a_(6,4)=1859/4104;
 a_(6,5)=-11/40; 
 % 4th order b-coefficients
 b4_(1)=25/216; 
 b4_(2)=0; 
 b4_(3)=1408/2565; 
 b4_(4)=2197/4104;
 b4_(5)=-1/5;
 % 5th order b-coefficients
 b5_(1)=16/135; 
 b5_(2)=0; 
 b5_(3)=6656/12825; 
 b5_(4)=28561/56430;
 b5_(5)=-9/50; 
 b5_(6)=2/55;
 for i=1:6
  c_(i)=sum(a_(i,:));
 end

% Initialization
t0 = tspan(1);
tfinal = tspan(2);
t = t0;
hmax = (tfinal - t)/2.5;
hmin = (tfinal - t)/1e9;
h = (tfinal - t)/100; % initial guess at a step size
y = y0(:);            % this always creates a column vector, y
tout = t;             % first output time
yout = y.';           % first output solution

if count==1,
 global counter
 if ~exist('counter'),counter=0; end
end % if count

if trace
 clc, t, h, y
end

% The main loop
   while (t < tfinal) & (h >= hmin)
      if t + h > tfinal, h = tfinal - t; end

      % compute the slopes
      k_(:,1)=feval(F,t,y);
      k_(:,2)=feval(F,t+c_(2)*h,y+h*(a_(2,1)*k_(:,1)));
      k_(:,3)=feval(F,t+c_(3)*h,y+h*(a_(3,1)*k_(:,1)+a_(3,2)*k_(:,2)));
     
      k_(:,4)=feval(F,t+c_(4)*h,y+h*(a_(4,1)*k_(:,1)+a_(4,2)*k_(:,2)+a_(4,3)*k_(:,3)));
      k_(:,5)=feval(F,t+c_(5)*h,y+h*(a_(5,1)*k_(:,1)+a_(5,2)*k_(:,2)+a_(5,3)*k_(:,3)+a_(5,4)*k_(:,4)));
      k_(:,6)=feval(F,t+c_(6)*h,y+h*(a_(6,1)*k_(:,1)+a_(6,2)*k_(:,2)+a_(6,3)*k_(:,3)+a_(6,4)*k_(:,4)+a_(6,5)*k_(:,5)));

      % increment the counter
      if count==1,
       counter = counter + 6;
      end % if

      % compute the 4th order estimate
      y4 = y + h*(b4_(1)*k_(:,1) + b4_(3)*k_(:,3) + b4_(4)*k_(:,4) + b4_(5)*k_(:,5));
      % compute the 5th order estimate
      y5 = y + h*(b5_(1)*k_(:,1) + b5_(3)*k_(:,3) + b5_(4)*k_(:,4) + b5_(5)*k_(:,5) + b5_(6)*k_(:,6));

      % estimate the local truncation error
      gamma1 = y5 - y4;

      % Estimate the error and the acceptable error
      delta = norm(gamma1,'inf');       % actual error
      tau = tol*max(norm(y,'inf'),1.0); % allowable error

      % Update the solution only if the error is acceptable
      if delta <= tau
         t = t + h;
         y = y5;
         tout = [tout; t];
         yout = [yout; y.'];
      end
      if trace
         home, t, h, y
      end

      % Update the step size
      if delta ~= 0.0
         h = min(hmax, 0.8*h*(tau/delta)^pow);
      end
   end;

   if (t < tfinal)
      disp('Step size grew too small.')
      t, h, y
   end
   return

