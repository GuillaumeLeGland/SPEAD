function [ADV] = jamstecrest_SinkingAdvection(C0,deltaz,wspeed,ndepths)

%========================================================================
%........................................................................
C0 = C0(:); %make it column vector (0 - 200m)
%........................................................................
CC0 = [0;C0;0]; %Add two boundary conditions with zero concentration (surface and bottom)
%........................................................................
%========================================================================
%........................................................................
J = [1:ndepths]; %depth indices
%........................................................................
JJ = J + 1; %Modify depth indices to match the modified boundary condtion [ie. CCO(JJ) = CO(J) ]
%........................................................................
%Checking test:
conc0 = C0(J);
conc0bis = CC0(JJ);
if conc0 ~= conc0bis
    'C0(J) is:'
    C0(J)
    'CC0(JJ) is:'
    CC0(JJ)
    error('C0(J) and CC0(JJ) should be the same!')
end
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FLUXES TO/FROM ALL GRID NODES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
GoesInn = wspeed(:).*CC0(JJ-1);
GoesOut = wspeed(:).*CC0(JJ);
%........................................................................
%GoesOut(end) = 0; %Force the bottom to be a "reflective" boundary condition (i.e. like hard rock)
% Sinking below 200m now allowed, balance by DIN input (Le Gland, 11/12/2019)
%........................................................................
%========================================================================

%%%%%%%%%%%%%%%%%%%
%ADVECTIVE SINKING:
%%%%%%%%%%%%%%%%%%%
ADV = (GoesInn - GoesOut) ./ deltaz; %First order up-wind method for advective sinking.

%%%%%%%%
%OUTPUT:
%%%%%%%%
ADV = ADV(:); %make it column vector.
