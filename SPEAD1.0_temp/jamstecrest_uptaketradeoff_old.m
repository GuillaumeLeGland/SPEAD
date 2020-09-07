function [mup,knp] = jamstecrest_uptaketradeoff(xaxis,keyChooseTradeoff)
global alp0 mup0 knp0 
global aalp amup aknp 

%===================================================================================
%...................................................................................
ESD = exp(xaxis); 
ESDmin = min(ESD);
ESDmax = max(ESD);
%...................................................................................
%===================================================================================
%MY APPROACH:
%...................................................................................
mup1 = 2.0; %[d-1]
knp1 = 0.5; %[mmolN*m-3]
alp1 = mup1/knp1; %[m3*mmolN-1*d-1]
Cte  = mup1*alp1; %[m3*mmolN-1*d-2]
%...................................................................................
xmax = max(xaxis(:));
%%xmax = log(200.0); %log([um])
betax = 0.02; 
%%betax = 0.05; 
%%betax = 4 / max(exp(xaxis(:))); %betax = 0.02 for xmax = log(200.0)
%...................................................................................
% $$$ mupmin = 0.2;
% $$$ mupmax = 4.0;
% $$$ mup001 = mupmax - (mupmax - mupmin)*(ESD-ESDmin) ./ (ESDmax-ESDmin); %Enforce mupmax y mupmin values.
%...................................................................................
% $$$ mupmax = 10.0;
% $$$ mupstar = (ESD-ESDmin) ./ (ESDmax-ESDmin); %Linear variable between [0 - 1] [n.d.] 
% $$$ mup001 = mupstar*mupmax;
%...................................................................................
% $$$ mup001 = mupmax * (1 - exp(-betax*ESD)); %Large cells grow faster.
% $$$ %%mup001 = mupmax * exp(-betax*ESD); %Small cells grow faster. 
% $$$ %%mup001 = mupmax * exp(-0.5*(xaxis-min(xaxis))); %Small cells faster.
mup001 = mup0 * exp(amup*xaxis);
%...................................................................................
knp001 = (mup001.^2)/Cte;
%...................................................................................
%===================================================================================
%LAN SMITH APPROACH:
%------------------------------------------------------------------------
%NOTE: To keep the same DIN trade-off, both "mup" and "alp" ** must ** be 
%multiplied by the environmental limitation factor (eg. Qpar or Qsst) 
%Otherwise, if Qpar or Qsst only multiplies "mup", the optimal size for a 
%given DIN value will shift up and down instead of remaining always at 
%the same ESDphy value. Multiplying Qpar or Qsst by both mup and alp
%means that kdin = mup/alp reamains invariant (ie. does not change with
%Qpar or Qss). Do *NEVER* multiply "kdin" by Qsst or Qpar!!!
%------------------------------------------------------------------------
%NOTE: To have a growth-affinity trade-off (i.e. mup x affinity = Cte)
%both "amup" and "aalp" ** must ** have the same value but with opposite
%sign ("amup" is positive and "aalp" is negative). 
%------------------------------------------------------------------------
%...................................................................................
factor = 1.0; %(e.g. full light conditions)
%%factor = 0.1; %(e.g. lower alfa and mup due to light limitation)
%...................................................................................
alp0factor = alp0*factor;
mup0factor = mup0*factor;
%...................................................................................
knp0factor = (mup0factor/alp0factor); %(unnafected by the multiplication factor)
%...................................................................................
% $$$ mup0factor = mup0*0.1; %(e.g. lower mup due to light limitation)
% $$$ knp0factor = knp0*1.9; %(e.g. higher knp0 due to light limitation)
%...................................................................................
alp002 = alp0factor * exp(aalp*xaxis);
mup002 = mup0factor * exp(amup*xaxis);
knp002 = knp0factor * exp(aknp*xaxis);
%...................................................................................
% $$$ figure(1000)
% $$$ subplot(2,2,1)
% $$$ plot(xaxis,alp002)
% $$$ title('alp')
% $$$ grid on
% $$$ subplot(2,2,2)
% $$$ plot(xaxis,mup002)
% $$$ title('mup')
% $$$ grid on
% $$$ subplot(2,2,3)
% $$$ plot(xaxis,knp002)
% $$$ title('knp')
% $$$ grid on
% $$$ return 
%...................................................................................
%===================================================================================
%OUTPUT:
%...................................................................................
if strcmp(keyChooseTradeoff,'Vallina')
    mup = mup001;
    knp = knp001;
elseif strcmp(keyChooseTradeoff,'Lanimal')
    mup = mup002;
    knp = knp002;
end
%***********************************************************************************
return
%===================================================================================
%PLOTS:
%...................................................................................
alp001 = mup001./knp001; %[m3*mmolN-1*d-1]
alp002 = mup002./knp002; %[m3*mmolN-1*d-1]
%...................................................................................
mupalp001 = mup001 .* alp001; 
mupalp002 = mup002 .* alp002; 
%...................................................................................
figure(1)
subplot(2,2,1)
plot(ESD,mup001)
xlabel('cell size')
ylabel('max growth')
title('mup')
grid on
subplot(2,2,2)
plot(ESD,knp001)
xlabel('cell size')
ylabel('half-sat')
title('knp')
grid on
%...................................................................................
subplot(2,2,3)
plot(ESD,mup002)
xlabel('cell size')
ylabel('max growth')
title('mup')
grid on
subplot(2,2,4)
plot(ESD,knp002)
xlabel('cell size')
ylabel('half-sat')
title('knp')
grid on
%...................................................................................
figure(2)
subplot(2,2,1)
plot(ESD,alp001)
xlabel('cell size')
xlabel('alp')
title('affinity')
grid on
subplot(2,2,2)
plot(ESD,alp002)
xlabel('cell size')
xlabel('alp')
title('affinity')
subplot(2,2,3)
plot(ESD,mupalp001)
xlabel('cell size')
ylabel('mup x alp')
title('growth x affinity')
subplot(2,2,4)
plot(ESD,mupalp002)
xlabel('cell size')
ylabel('mup x alp')
title('growth x affinity')
%...................................................................................
%===================================================================================
%***********************************************************************************
return

%===================================================================================
%...................................................................................
mup1 = 2.0; %[d-1]
knp1 = 0.5; %[mmolN*m-3]
alp1 = mup1/knp1; %[m3*mmolN-1*d-1]
Cte  = mup1*alp1; %[m3*mmolN-1*d-2]
%...................................................................................
MUP001 = [1.25, 1.75, 2.25, 2.75];
KNP001 = (MUP001.^2)/Cte; 
%...................................................................................
KNP002 = [0.15, 0.30, 0.60, 1.20];
MUP002 = sqrt(Cte*KNP002); 
%...................................................................................
%===================================================================================
%...................................................................................
mup0 = 1.25; 
pcnt = (1/3);
MUP = [mup0];
for i = 1:3 
    mupi = mup0 + (pcnt*mup0) 
    MUP = [MUP,mupi]
    mup0 = mupi;
end
MUP003 = MUP;
KNP003 = (MUP003.^2)/Cte; 
%...................................................................................
%===================================================================================
