function [ynoise]=myStochasticNoise(ndays,pcntmax)
    
%..............................
time=[1:ndays];
%..............................

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GENERATE STOCHASTIC NOISE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%..............................
% $$$ dtnoise=1;
dtnoise=5;
%..............................
tnoise=[dtnoise:dtnoise:ndays];
mnoise=length(tnoise);
%..............................
xnoise=(1-pcntmax)+(rand(1,mnoise)*(2*pcntmax));
%..............................
% $$$ plot(time,xnoise,'b-',time,xnoise,'r*')
%..............................

% $$$ return

%%%%%%%%%%%%%%%%%%%%
%INTERP TIME SERIES:
%%%%%%%%%%%%%%%%%%%%
xnoiseI=interp1(tnoise,xnoise,time);

%%%%%%%%%%%%%%%%%%%%
%SMOOTH TIME SERIES:
%%%%%%%%%%%%%%%%%%%%
win=3;
[wnoise,junk]=myrunmean(xnoiseI,win,'Temporal','Periodica');
%..............................
% $$$ plot(wnoise,'*-')
%..............................

%%%%%%%%
%OUTPUT:
%%%%%%%%
%.............
% $$$ ynoise=xnoise;
%.............
ynoise=wnoise;
%.............

%%%%%%
%PLOT:
%%%%%%
% $$$ plot(ynoise,'*-')
