function [UML,UMLnoisy]=mypopovaUML(it1,it2,hmin1,hmax1)

T=365;
tspan=[1:T];

%CALCULATE UML:
for i=1:it1        %(ie. t=1 - t=118)
    h(i)=-(hmax1-hmin1)*(sin(pi*(i+T)/T)-1) / (sin(pi*it1/T)+1)+hmin1;
end

for i=it2:365      %(ie. t=178 - t=365)
    h(i)=-(hmax1-hmin1)*(sin(pi*i/T)-1) / (sin(pi*it1/T)+1)+hmin1;
end

for i=it1+1:it2-1 %(ie. t=119 - t=177) 
    h(i)=hmax1-(i-it1)*(h(it1)-h(it2))/(it2-it1); %usa las h(it1) y h(it2) obtenidas en los otros bucles DO
end
UML=h;

%ADD NOISE:
%a) Generate random variability de hasta el 50%:
randvar=rand(1,365)/2; 

%b) Generate whether it is positive variability (MLD increses) or neg.var.:
posnegmask=ones(1,365);
xrandom=rand(1,365);
I=find(xrandom<0.50);
posnegmask(I)=-1;
randmask=randvar.*posnegmask;
UMLnoisy=UML+(UML.*randmask);

%..............................
% $$$ figure(1)
% $$$ plot(tspan,UML*(-1),'-r.',tspan,UMLnoisy*(-1),'-b.')
% $$$ legend('UML1','UML2')
% $$$ axis([1 365, -300 0])
%..............................


