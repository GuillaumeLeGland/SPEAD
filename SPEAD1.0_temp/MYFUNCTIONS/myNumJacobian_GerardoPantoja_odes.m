function [Xdot] = myNumJacobian_GerardoPantoja_odes(ti,X0)
global a b c d 
global mp mz

N = X0(1);
P = X0(2);
Z = X0(3);

Ndot =   -a*P*N + (1-c)*b*P*Z + d*Z + mp*P^2 + mz*Z^2;
Pdot =    a*P*N -       b*P*Z       - mp*P^2;
Zdot =               c *b*P*Z - d*Z - mz*Z^2;
    
Xdot(1) = Ndot;
Xdot(2) = Pdot;
Xdot(3) = Zdot;
Xdot = Xdot(:);

return

