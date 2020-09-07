function [retval] = isOctave()
%*******************************************************************
% Use: [retval] = isOctave()
% Returns: [true] if the environment is Octave.
% Returns: [retval = 1] (trues) if the operating environment is Octave.
% Returns: [retval = 0] (false) if the operating environment is Matlab.
%*******************************************************************
%===================================================================
%...................................................................
  persistent cacheval;  % speeds up repeated calls
  if isempty (cacheval)
    cacheval = (exist ('OCTAVE_VERSION', 'builtin') == 5);
  end
  retval = cacheval;
%...................................................................
% <https://www.mathworks.com/matlabcentral/fileexchange/28847-is-octave>
% If OCTAVE_VERSION is a built-in function, then we must be in Octave.
% Since the result cannot change between function calls, it is cached in a
% persistent variable. isOctave cannot be a persistent variable, because it
% is the return value of the function, so instead the persistent result must
% be cached in a separate variable.
%===================================================================
%*******************************************************************
return
