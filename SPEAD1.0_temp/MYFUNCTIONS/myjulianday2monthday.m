function [month,mday]=myjulianday2monthday(Jday)
%**********************************
%Syntaxis: [month,mday]=myjulianday2monthday(jday)
%**********************************
%................................
% $$$ jday=224; %(12 august)
% $$$ %jday=31;
%................................

n=length(Jday);
x=[0,31,59,90,120,151,181,212,243,273,304,334,365];
month=[];
mday=[];
for i=1:n
    Jdayi=Jday(i);
    I=find(Jdayi<=x);
    monthi=I(1)-1; %month para el dia juliano i.
    mdayi=Jdayi-x(monthi); %day of the month para el dia juliano i.
    %STOCKAGE:
    month=[month;monthi];
    mday=[mday;mdayi];
end

% $$$ [jday,month,mday]
 
 
