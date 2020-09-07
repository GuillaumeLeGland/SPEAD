function [MLD,MLDnoisy]=mypopovaMLD(ndays,it1,it2,hmin1,hmax1)

T=ndays;
tspan=[1:T];

%CALCULATE MLD:
for i=1:it1        %(ie. t=1 - t=118)
    h(i)=-(hmax1-hmin1)*(sin(pi*(i+T)/T)-1) / (sin(pi*it1/T)+1)+hmin1;
end

for i=it2:ndays      %(ie. t=178 - t=365)
    h(i)=-(hmax1-hmin1)*(sin(pi*i/T)-1) / (sin(pi*it1/T)+1)+hmin1;
end

for i=it1+1:it2-1 %(ie. t=119 - t=177) 
    h(i)=hmax1-(i-it1)*(h(it1)-h(it2))/(it2-it1); %usa las h(it1) y h(it2) obtenidas en los otros bucles DO
end
MLD=h;

%ADD NOISE:
%a) Generate random variability de hasta el 50%:
randvar=rand(1,ndays)/2; 

%b) Generate whether it is positive variability (MLD increses) or neg.var.:
posnegmask=ones(1,ndays);
xrandom=rand(1,ndays);
I=find(xrandom<0.50);
posnegmask(I)=-1;
randmask=randvar.*posnegmask;
MLDnoisy=MLD+(MLD.*randmask);

%..............................
% $$$ figure(1)
% $$$ plot(tspan,MLD*(-1),'-r.',tspan,MLDnoisy*(-1),'-b.')
% $$$ legend('MLD1','MLD2')
% $$$ axis([1 ndays, -300 0])
%..............................


