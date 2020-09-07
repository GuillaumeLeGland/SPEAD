function [EDAYuml]=mySRD(EDAY0,MLD,KEXT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INTEGRO LA EDAY EN LA PROFUNDIDAD PARA CADA PIXEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------
%NOTA:
%a) INTEGRAL de PAR+ EN Z(m) (sol.analitica): 
%integral de "I(z) = Iz = Io*exp(-a*z)" respecto de z(m) = int(Iz)dz.
%SOL: int(Iz)dz = (Io/a)*(1-exp(-a*z))
%---------------------------------------------
IEDAYdz=(EDAY0./KEXT).*(1-exp(-KEXT.*MLD)); %int(Iz)dz; [W*m-2] x [m] = [W*m-1]
EDAYuml=IEDAYdz./MLD; %[W*m-1] x [m-1] = [W*m-2]
