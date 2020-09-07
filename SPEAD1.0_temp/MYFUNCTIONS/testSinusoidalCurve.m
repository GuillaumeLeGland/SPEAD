close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%
%GET SINUSOIDAL CURVES:
%%%%%%%%%%%%%%%%%%%%%%%
%............
T=365;
years=1;
tmax=T*years;
%............
ymin=20;
ymax=30;
%............
Ysin=[];
%............
ncurves=5;
%............
for i=1:ncurves
    %=================
    %SINUSOIDAL CURVE:
    %=================
    tau=(i/ncurves)*T;
    [ysin]=mySinusoidalCurve(T,tau,tmax,ymin,ymax);
    
    %=========
    %STOCKAGE:
    %=========
    Ysin=[Ysin;ysin];
end
sumYsin=sum(Ysin);

%%%%%%%%%%%
%ADD NOISE:
%%%%%%%%%%%
YsinNoisy=[];
for i=1:ncurves
    %=================
    %SINUSOIDAL CURVE:
    %=================
    ysin=Ysin(i,:);

    %=========
    %ADD NOISE:
    %=========
    [ysinNoisy]=myAddNoise(ysin);
    
    %=========
    %STOCKAGE:
    %=========
    YsinNoisy=[YsinNoisy;ysinNoisy];
end
sumYsinNoisy=sum(YsinNoisy);

%%%%%%%
%PLOTS:
%%%%%%%
figure(1)
subplot(2,2,1)
plot([1:tmax],Ysin)
legend('sin1','sin2','sin3','sin4','sin5')
subplot(2,2,3)
plot([1:tmax],sumYsin)
subplot(2,2,2)
plot([1:tmax],YsinNoisy)
legend('sin1','sin2','sin3','sin4','sin5')
subplot(2,2,4)
plot([1:tmax],sumYsinNoisy)

