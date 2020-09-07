function [TURNOVER,I,E,N] = mySpeciesTurnoverRateIndex(Mdata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAKE DATA AS PRESENCE/ABSENCE ONLY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================
%...................................................................
[nspecies,ntimes] = size(Mdata);
%...................................................................
Mzerones = zeros(nspecies,ntimes);
%...................................................................
Iones = find(Mdata > 0);
%...................................................................
Mzerones(Iones) = 1.0;
%...................................................................
M = Mzerones; 
%...................................................................
%===================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE TURNOVER RATE INDEX:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================
%-------------------------------------------------------------------
%A: sample-1 
%B: sample-2 
%-------------------------------------------------------------------
%TURNOVER = (I + E) / N 
%
% where:
%
% I :: Invasions.
% E :: Extinctions.
% N :: Number of species.
%-------------------------------------------------------------------
%===================================================================
for iday = 1:ntimes
    A = M(:,iday);
    for jday = 1:ntimes
	B = M(:,jday);

	N01 = length(find(A==0 & B==1)); %Invasions.
	N10 = length(find(A==1 & B==0)); %Extinctions.
% $$$ 	N11 = length(find(A==1 | B==1)); %Species number. %WRONG!!!!

	A1 = length(find(A==1)); %Species number.
	B1 = length(find(B==1)); %Species number.
	N11 = A1 + B1; %Okay.

	A1bis = sum(A==1); %Species number.
	B1bis = sum(B==1); %Species number.
	N11bis = A1bis + B1bis; %Okay.

% $$$ 	A1,A1bis
% $$$ 	B1,B1bis
% $$$ 	N11,N11bis,N11tris
% $$$ 	pause 

	I(iday,jday) = N01; %Invasions.
	E(iday,jday) = N10; %Extinctions.
	N(iday,jday) = N11; %Species number.

    end
end
%...................................................................
TURNOVER = (I + E) ./ N;
%...................................................................
%===================================================================
%*******************************************************************
return
