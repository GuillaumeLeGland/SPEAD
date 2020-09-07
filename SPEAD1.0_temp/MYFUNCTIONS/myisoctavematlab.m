function [uiIsOctave,uiIsMatLab] = myisoctavematlab
% <https://www.mathworks.com/matlabcentral/fileexchange/23868-is-this-matlab-or-octave>
% if uimatlab 
% disp('This is MATLAB.') 
% elseif uioctave 
% disp('This is Octave.') 
% else 
% disp('I do not know what user environment this is.') 
% end

uiIsOctave = false;
uiIsMatLab = false;
LIC = license('inuse');
for elem = 1:numel(LIC)
    envStr = LIC(elem).feature;
    if     strcmpi(envStr,'octave')
        uiIsOctave = true;
    elseif strcmpi(envStr,'matlab')
        uiIsMatLab = true;
    end
end
return

% ----
% <https://www.mathworks.com/matlabcentral/fileexchange/23868-is-this-matlab-or-octave>
% function uiIsOctave = uioctave
% uiIsOctave = false;
% LIC = license('inuse');
% for elem = 1:numel(LIC)
%     envStr = LIC(elem).feature;
%     if strcmpi(envStr,'octave')
%         uiIsOctave = true;
%         break
%     end
% end
% ----
% MATLAB OCTAVE CHECK:
% LIC = license('inuse');
% envStr = LIC(1).feature;
% if     strcmpi(envStr,'matlab')
%     addpath(genpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYTOOLBOX/'))
%     opengl('neverselect')
% elseif strcmpi(envStr,'octave')
%     % <https://octave.sourceforge.io/symbolic/function/@sym/pretty.html>
%     pkg('load','statistics');
%     pkg('load','symbolic')
%     pkg('load','odepkg')
%     sympref('display','flat') %USE THIS ONE BETTER
%     %%sympref('display','ascii')
% end
% ----
