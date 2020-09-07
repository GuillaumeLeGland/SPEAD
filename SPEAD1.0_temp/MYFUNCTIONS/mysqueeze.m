function [Aout] = mysqueeze(A,squeezedim)

%************************************************************************
%------------------------------------------------------------------------
%NOTE: Squeeze *only* along the SELECTED dimension (not along all of them)
%This function "mysqueeze.m" modifies the original "squeeze.m" in order
%to force squeezing along a *selected* singleton dimension exclusively 
%instead that along all the singleton dimensions. 
%------------------------------------------------------------------------
%Use: [Aout] = mysqueeze(A,squeezedim)
% 
% # A = Array (2D or larger)
% # Aout = Array squeezed along a selected dimension exclusively. 
% # squeezedim = Selected dimension to perform the squeeze alone. 
%************************************************************************
%........................................................................ 
Asize = size(A);
Zerones = zeros(1,length(Asize));
Zerones(squeezedim) = 1; %Tag 'squeezedim' as "one" (otherwise are "zeros")
U = (Zerones==1); %Pointer to 'squeezedim' exclusively.
Asize(U) = []; % Remove singleton dimension for 'squeezedim' exclusively.
Asize = [Asize, ones(1,2-length(Asize))]; % Make sure Asize is at least 2-D
Aout = reshape(A,Asize);
%........................................................................ 

%************************************************************************
return
%........................................................................ 
Asize = size(A),pause
Asize(Asize==1) = []; % Remove singleton dimensions.
Asize = [Asize, ones(1,2-length(Asize))]; % Make sure siz is at least 2-D
AoutBis = reshape(A,Asize);
%........................................................................ 
% <http://stackoverflow.com/questions/23703419/squeeze-some-of-singleton-dimensions-in-matlab> 
Asize = size(A);
U = (Asize ~= 1); 
U(squeezedim) = 1;
Aout = reshape(A,Asize(U)); 
%........................................................................ 


