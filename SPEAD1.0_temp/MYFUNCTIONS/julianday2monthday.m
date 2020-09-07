function [month,mday]=julianday2monthday(jday)
%**********************************
%Syntaxis: [month,mday]=julianday2monthday(jday)
%**********************************
% $$$ jday=224; %(12 august)
% $$$ %jday=31;

x=[0,31,59,90,120,151,181,212,243,273,304,334,365];
I=find(jday<=x);
month=I(1)-1;
mday=jday-x(month);

% $$$ [jday,month,mday]
 
 
