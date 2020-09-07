function [mX,tspan]=myDailyAverages(X,dt,varargin)

%*************************************
%Use: [mX,tspan]=myDailyAverages(X,dt)
%*************************************

%.....................    
windowsize='Daily';
%.....................
nvarargin=length(varargin);
if nvarargin==1
    windowsize=varargin{1};
end
%.....................
[m,n]=size(X);
if strcmp(windowsize,'Daily')
    daywindow=1/dt; %One point per day.
elseif strcmp(windowsize,'Subdaily')
    daywindow=0.25/dt; %Four points per day.
end
%.....................
% $$$ windowsize
% $$$ daywindow %It has to be an integer.
%.....................
%Check:
xmod = mod(daywindow,1); %must be zero.
if xmod ~= 0
    xmod
    disp('Error! "dt" must be an integer multiple of 1!! (eg. 0.05, 0.1, etc)')
    pause
end
%.....................
mX=[];
stdX=[];
ispan=1:daywindow:n;
%%tspan=ispan/daywindow;
%.....................
for i=ispan
    Xi=X(:,i:i+(daywindow-1));
    
    mXi=mean(Xi,2);
    stdXi=std(Xi,0,2);

    mX=[mX,mXi]; %(9vars,365days)
    stdX=[stdX,stdXi];
end
%.....................

%TIME OUTPUT:
t0=1*dt; %t0=dt
tmax=n*dt; %tmax=365
dtspan=daywindow*dt;
tspan=[t0:dtspan:tmax];

