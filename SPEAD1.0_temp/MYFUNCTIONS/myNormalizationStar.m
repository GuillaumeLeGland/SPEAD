function [Astar]=myNormalizationStar(A,varargin)

nvarargin = length(varargin);
if nvarargin == 1
    Aminmax = varargin{1};
    Amin = Aminmax(1);
    Amax = Aminmax(2);
else
    Amin = min(A(:));
    Amax = max(A(:));
end

Astar = (A - Amin) / (Amax - Amin); %[0 - 1]
