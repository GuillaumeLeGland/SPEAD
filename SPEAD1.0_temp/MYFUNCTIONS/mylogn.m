function [y] = mylogn(x,Base)
%************************************************************************
% <http://www.mathworks.com/matlabcentral/newsreader/view_thread/172539/>
% 
% Use: [y] = mylogn(x,Base);
%
%************************************************************************

y = log10(x)/log10(Base); 

return

%************************************************************************
% <http://www.rapidtables.com/math/algebra/logarithm/Logarithm_Base_Change.htm> 
%------------------------------------------------------------------------
%NOTE:
% logb(x) = loga(x)  / loga(b);
% log(x)  = log10(x) / log10(exp(1));
%------------------------------------------------------------------------

