function [ESDphyAveUpper,ESDphyAveLower,ESDphyAve,ESDphyStd,ESDphyCV] = jamstecrest_meansigmaESD(logESDphyAve,logESDphyStd)

%===================================================================================
%-----------------------------------------------------------------------------------
%NOTE: 
%Doing ESDave = exp(logESDave) is okay. 
%Doing ESDstd = exp(logESDstd) is wrong!!!! (see how to do it right below) 
%-----------------------------------------------------------------------------------
%...................................................................................
%COMPUTING AVERAGE AND RANGE OF CELL SIZE IN ABSOLUTE SCALE FROM LOG(CELL SIZE) SCALE:
ESDphyAve = exp(logESDphyAve); %Okay (i.e. ESDave = median(ESD))
ESDphyAveUpper = ESDphyAve .* exp(logESDphyStd);
ESDphyAveLower = ESDphyAve ./ exp(logESDphyStd);
%...................................................................................
%COMPUTING COEFFICIENT OF VARIATION OF CELL SIZE IN ABSOLUTE SCALE FROM LOG(CELL SIZE) SCALE:
%(ERROR! ALL THIS SHIT IS *NOT* WORKING!!!!! I DON'T KNOW WHY... CHECK IT OUT)
ESDphyStd = (ESDphyAveUpper - ESDphyAveLower); %Crude estimation but should be very accurate.
ESDphyStdChap = (ESDphyAve ./ exp(logESDphyStd)) .* (exp(logESDphyStd).^2 - 1d0); %Theoretical equation. 
ESDphyCV = ESDphyStd./ESDphyAve; %Coefficient of variation (std/mean)
%...................................................................................
%===================================================================================
