close all
clear all

%f(1) =  a*P - b*P*Z - mp*P^2;
%f(2) = -b*P*Z - c*Z - mz*Z^2;

global a b c d 
global mp mz

myfun = '@myNumJacobian_GerardoPantoja_odes';

mp = 0.1;
mz = 0.2;
a = 1.2;
b = 0.3;
c = 0.4;
d = 0.05;

N = 0.50;
P = 0.25;
Z = 0.25;

x0(1) = N;
x0(2) = P;
x0(3) = Z;
x0 = x0(:);

Ntot0 = N + P + Z;

dt1 = 1/24;
tmax = 90;
tspan1 = [dt1:dt1:tmax]

%EULER:
Nout = [];
Pout = [];
Zout = [];
for ti = tspan1(:)'

    Ndot =   -a*P*N + (1-c)*b*P*Z + d*Z + mp*P^2 + mz*Z^2;
    Pdot =    a*P*N -       b*P*Z       - mp*P^2;
    Zdot =               c *b*P*Z - d*Z - mz*Z^2;

    N = N + (Ndot*dt1);
    P = P + (Pdot*dt1);
    Z = Z + (Zdot*dt1);

    Nout = [Nout,N];
    Pout = [Pout,P];
    Zout = [Zout,Z];
end
Ntot_out = Pout + Zout + Nout;

Nssp = Nout(end);
Pssp = Pout(end);
Zssp = Zout(end);

%RUNGE KUTTA:
dt2 = 1;
tmax = 90;
tspan2 = [dt2:dt2:tmax]

[tode,Xode] = ode45(eval(myfun),tspan2,x0); %RK04 with adaptative time-step solution.
Xode = Xode';
Node = Xode(1,:);
Pode = Xode(2,:);
Zode = Xode(3,:);
Ntot_ode = Pode + Zode + Node;

Nssp = Node(end);
Pssp = Pode(end);
Zssp = Zode(end);

Xssp = [Nssp,Pssp,Zssp];

%PLOTS:

figure(1)
subplot(2,2,1)
plot(tspan1,Nout,'r-')
hold on
plot(tspan1,Pout,'g-')
hold on
plot(tspan1,Zout,'k-')
hold on
plot(tspan1,Ntot_out,'b-')
hold off

subplot(2,2,2)
plot(tspan2,Node,'r-')
hold on
plot(tspan2,Pode,'g-')
hold on
plot(tspan2,Zode,'k-')
hold on
plot(tspan2,Ntot_ode,'b-')
hold off

%JACOBIAN:
[Jf] = myNumJacobian(eval(myfun),Xssp);
%%[Jf] = myNumJacobianOctave(eval(myfun),Xssp);
%%[Jf] = myNumJacobianFast(eval(myfun),Xssp)

X  = Xssp;
m  = length(X);
Jfbis = zeros(m,m);
dx = 0.001;
for j = 1:m
    %......................
    X1    = X;
    X2    = X;
    %......................
    X2(j) = X(j) + dx;
    %......................
    f1 = feval(eval(myfun),dt1,X1); %This way allways works ok.
    f2 = feval(eval(myfun),dt1,X2);
    %......................
    Jfbis(:,j) = (f2 - f1)/dx;  % partial derivatives in j-th row.
    %......................
end
eig(Jf)

return
