function [Ysmooth] = myConvRunMean(Y,nptos)

%****************************************
%[Ysmooth]=myConvRunMean(Y,nptos)
%
% Y: Data vector original.
% Ysmooth: Data vector smoothed.
% nptos: Running window length-size.
%****************************************
%============================================================
%............................................................
Y = Y(:); 
%............................................................
%%keySignalFrequency = 'nan';
keySignalFrequency = 'Periodica';
%............................................................
msize = length(Y); 
WindowSize = nptos; 
Win = ones(1,WindowSize)/WindowSize; 
%............................................................
%%shape = 'full'; %DONT USE!!!! SALE FATAL.
shape = 'same'; %OKAY 
%............................................................
if strcmp(keySignalFrequency,'Periodica')
    YE = [Y;Y;Y]; %Extended vector by adding before and after. 
    YEsmooth  = conv(YE,Win,shape); 
    Ysmooth   = YEsmooth(msize+1:end-msize);
else
    Ysmooth   = conv(Y,Win,shape); 
end 
%............................................................
nsize = length(Ysmooth);
%............................................................
% $$$ msize
% $$$ nsize
% $$$ pause 
%============================================================
return
%************************************************************
%------------------------------------------------------------
%<http://stackoverflow.com/questions/3453663/computing-running-averages-in-matlab>
%............................................................
ysin = sin([0:(pi/24):(2*pi)])
%............................................................
pcnt = 1.0;
msize = length(ysin);
xrandom = (0.5 - rand(1,msize))*pcnt;
%............................................................
ysinNoisy = ysin + (max(ysin).*xrandom);
%............................................................
windowSize = 5;
Win = ones(1,windowSize)/windowSize;
%............................................................
Y = ysinNoisy;
%............................................................
Ysmooth = conv(Y,Win); %If you dont use the 'valid' argument, Ysmooth is larger than Y. To correct:
%............................................................
halfSize = floor(windowSize/2);
Ysmooth = Ysmooth([(halfSize+1):(end-halfSize)]);
%............................................................
YsmoothBis  = conv(Y,Win,'same');
YsmoothTris = conv(Y,Win,'valid');
%............................................................
figure(10)
subplot(2,2,1)
plot([1:msize],ysin,'*r-',[1:msize],ysinNoisy,'.b-',[1:msize],Ysmooth,'+g-')
%............................................................
%------------------------------------------------------------
