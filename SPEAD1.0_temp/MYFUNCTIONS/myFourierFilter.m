function [SignalLowPass,SignalHighPass]=myFourierFilter(SignalNoisy)

%********************************************
%APLICO UN FILTRO DE TRANSFORMADA DE FOURIER:
%********************************************
%.....................
n=length(SignalNoisy);
if mod(n,2)==1 %n is a odd number.
    SignalNoisy=SignalNoisy(1:end-1);
end
n=length(SignalNoisy); %Has to be an even number!!
%.....................

%=======================================================================
%REPEAT 1st AND last DATA POINT TO AVOID PROBLES WITH NON-PERIDI SIGNAL:
%=======================================================================
nquarter=n/4;
x1=SignalNoisy(1);
xend=SignalNoisy(end);
X1=[];
Xend=[];
for i=1:nquarter
    X1=[X1,x1];
    Xend=[Xend,xend];
end
SSignalNoisy = [X1,SignalNoisy,Xend];
nn=length(SSignalNoisy);

%................
% $$$ plot(SSignalNoisy)
% $$$ pause
%................

%====================
%CONSTRUYO EL FILTRO:
%====================
%.....................
% $$$ coef=10;
% $$$ coef=20;
% $$$ coef=30; %USAR ESTE.
coef=40;%  inicial 40
%.....................
nhalf=nn/2;
%.....................
yhalf1=[1:nhalf];
yhalf2=[nhalf:-1:1];
Y=[yhalf1,yhalf2];
%.....................
LowPassFilter=1./(1+(Y/coef).^10); %paso bajo
HighPassFilter=1-LowPassFilter; %paso alto.
%.....................

%=====================================================
%OBTENGO LA TRANSFORMAD DE FOURIER (FRECUENCY DOMAIN):
%=====================================================
FFT=fft(SSignalNoisy); %Transformada de Fourier

%=============================================
%APLICO EL FILTRO A LA TRANSFORMAD DE FOURIER:
%=============================================
FFTlowpass=FFT.*LowPassFilter; %Filtro altas frecuencias.
FFThighpass=FFT.*HighPassFilter; %Filtro bajas frecuencias.

%====================
%BACK TO TIME DOMAIN:
%====================
SSignalLowPass=real(ifft(FFTlowpass)); 
SSignalHighPass=real(ifft(FFThighpass)); 

%%%%%%%%%%
%GRAFICAS:
%%%%%%%%%%
figure(1)
plot(SSignalNoisy,'r')
hold on
plot(SSignalLowPass,'b-') %High Frecuency Filtered.
hold on
plot(SSignalHighPass,'g-') %Low Frecuency Filtered.

%%%%%%%%%
%OUTPUTS:
%%%%%%%%%
i1=nquarter+1;
i2=nquarter+n;
SignalLowPass=SSignalLowPass(i1:i2); %High Frecuency Filtered.
SignalHighPass=SSignalHighPass(i1:i2); %Low Frecuency Filtered.
