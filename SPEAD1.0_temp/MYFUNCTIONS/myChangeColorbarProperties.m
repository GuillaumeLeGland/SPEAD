function [newbarpos]=myChangeColorbarProperties(hc,gamma1,gamma2,gamma3,gamma4)

%**************************    
%This program changes the properties of the colorbar. The parameters
%gamma will multiply the original positions of the colorbar:
%
%newbarpos(1)=barpos(1)*gamma1; %left.
%newbarpos(2)=barpos(2)*gamma2; %bottom.
%newbarpos(3)=barpos(3)*gamma3; %width.
%newbarpos(4)=barpos(4)*gamma4; %height.
%
%and the new colorbar position will be used.
%**************************    

barpos=get(hc,'position');
newbarpos(1)=barpos(1)*gamma1; %left.
newbarpos(2)=barpos(2)*gamma2; %bottom.
newbarpos(3)=barpos(3)*gamma3; %width.
newbarpos(4)=barpos(4)*gamma4; %height.

%.......
I=find(newbarpos<=0); %in case some value is negative or zero.
newbarpos(I)=eps;
%.......

set(hc,'position',newbarpos);
