function [PARoutput]=myPARconversion(PARinput,type)
    
%***************************************************
%Use: [PARoutput]=myPARconversion(PARinput,type)
%
%type = 'Einstains-to-Watts'
%type = 'Watts-to-Einstains'
% 
% <http://www.egc.com/useful_info_lighting.php> 
%***************************************************

%.......................................................
% $$$ fcPAR0 = 1/4.6;         %conversion from [uEin*m-2*s-1] to [W*m-2] 
% $$$ fcPAR1 = 1/0.2174;      %conversion from [W*m-2] to [uEin*m-2*s-1] = 4.6 approx.
% $$$ fcPAR2 = 3600*24/1e6;   %conversion from [uEin*m-2*s-1] to [Ein*m-2*d-1].
% $$$ fcPAR3 = fcPAR1*fcPAR2; %conversion from [W*m-2] to [Ein*m-2*d-1].
% $$$ fcPARconv = 1/fcPAR3;   %conversion factor from [Einstein*m-2*d-1] to [W*m-2] at 400-700nm (PAR+)
%.......................................................
fcPAR = 2.5; %factor de conversion de [Einstein*m-2*d-1] a [W*m-2] para 400-700nm (PAR)
%.......................................................
if strcmp(type,'Einstains-to-Watts')
    PARoutput =   (fcPAR)*PARinput; %[W*m-2]
elseif strcmp(type,'Watts-to-Einstains')
    PARoutput = (1/fcPAR)*PARinput; %[Einstein*m-2*d-1]
end
%.......................................................
