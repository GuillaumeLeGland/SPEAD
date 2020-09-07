function [Wnoisy]=myAddNoise(W,varargin) %(see also "mypopovaUML.m")

m=length(W);    
    
%ADD NOISE:
%=========================================================
%1) Generate random variability de hasta el 50% o el 100%:
%=========================================================
%.........................
nvarargin=length(varargin);
%.........................
if nvarargin==1
    pcntmax=varargin{1}
else
% $$$     pcntmax=0.05;
% $$$     pcntmax=0.1;
    pcntmax=0.5; %USAR ESTE!!!
% $$$     pcntmax=1.0;
end
%.........................
randvar=rand(1,m)*pcntmax; 
%.........................

%==========================================================================
%2) Generate whether it is positive variability (W increses) or neg.var.:
%==========================================================================
posnegmask=ones(1,m);
%...............
xrandom=rand(1,m);
I=find(xrandom<0.5);
posnegmask(I)=-1; %Sometimes (+) and sometimes (-)
%...............
% $$$ posnegmask(:)=-1; %Always (-)
%...............
% $$$ posnegmask(:)=+1; %Always (+)
%...............

%===============
%3) Apply noise:
%===============
randmask=randvar.*posnegmask;

%æææææææææææææææææææææææææææææææææææææææææææææææææææ
%a) %random (+/-) 50 percent of value at each point:
%æææææææææææææææææææææææææææææææææææææææææææææææææææ
%.....................
% $$$ Wnoisy=W+(W.*randmask); 
%.....................

%æææææææææææææææææææææææææææææææææææææææææææææææææææææ
%b) Random (+/-) 50 percent of time-series mean value:
%æææææææææææææææææææææææææææææææææææææææææææææææææææææ
%.....................
meanW=mean(W(:));
Wnoisy=W+(meanW.*randmask);
%.....................

%..............................
% $$$ figure(100)
% $$$ plot([1:m],W,'-r.',[1:m],Wnoisy,'-b.')
% $$$ legend('W1','W2')
%..............................


