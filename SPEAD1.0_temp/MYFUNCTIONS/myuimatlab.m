function [uiIsMatLab] = uimatlab
% <https://www.mathworks.com/matlabcentral/fileexchange/23868-is-this-matlab-or-octave->
% if uimatlab 
% disp('This is MATLAB.') 
% elseif uioctave 
% disp('This is Octave.') 
% else 
% disp('I do not know what user environment this is.') 
% end

uiIsMatLab = false;
LIC = license('inuse');
for elem = 1:numel(LIC)
    envStr = LIC(elem).feature;
    if strcmpi(envStr,'matlab')
        uiIsMatLab = true;
        break
    end
end
return
