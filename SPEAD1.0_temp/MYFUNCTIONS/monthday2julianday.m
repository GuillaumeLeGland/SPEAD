function [Jday]=monthday2julianday(M,D)
%**********************************
%Syntaxis: [Jday]=julianday2monthday(M,D)
%**********************************
%M=mes;
%D=day;

%AUN NO FUNCIONA PARA VECTORES MATRICIALMENTE.
    
if min(size(D))>1
    error('M y D deben ser escalares o vectores')
end
m=length(D);    
if m==1 %escalar
    month=M;
    day=D;
%................................
% $$$ month=8;day=12; %jday=224;
% $$$ month=1;day=31; %jday=31;
%................................
x=[0,31,59,90,120,151,181,212,243,273,304,334,365];
jday=x(month)+day;
%................................
% $$$ [month,mday,jday]
%................................
Jday=jday;
else %vector
    %................................
% $$$     M=[1,8];D=[31,12]; %(31 Jan y 12 Ago)
    %................................
    %a) MATRICIALMENTE:
% $$$     xlim=[0,31,28,31,30,31,30,31,31,30,31,30,31];
% $$$     Xlim=ones(m,1)*xlim
% $$$     X=Xlim(:,1:M);
% $$$     Jday=sum(X,2)+D;
    %b) CON BUCLES:
    xlim=[0,31,28,31,30,31,30,31,31,30,31,30,31];
    Jday=[];
    for i=1:m
	i
	day=D(i);
	month=M(i);
	x=xlim(1:month);
	jday=sum(x)+day;
	Jday=[Jday,jday];
    end
end
Jday=Jday(:);















