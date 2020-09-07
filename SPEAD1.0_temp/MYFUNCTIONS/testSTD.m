clear all
close all

%You better see "egSTDforSeveralVariables.m" (it is a better and more complete example).

time=[1:365];

x1=rand(365,1); %Phyto.
x2=rand(365,1);
x3=rand(365,1);

y1=10+(rand(365,1)*2); %Zoo.
y2=10+(rand(365,1)*2);
y3=10+(rand(365,1)*2);

x=x1+x2+x3; %PhyTot.
y=y1+y2+y3; %ZooTot.

%%%%%%%%%%%%%
%INDIV. CONC:
%%%%%%%%%%%%%

%%%%%%%%%%%%%
%CONC. TOTAL:
%%%%%%%%%%%%%
%............
pvarx=var(x);
vary=var(y);
covxy=cov(x,y);
%............
sumCovXY=sum(covxy(:));
SqrtSumCovXY=sqrt(sumCovXY); %Output.
%............

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZE DATA, (X+Y) ADDITION, AFTER STD CALCULATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumxy=x+y;
meanSumXY=mean(sumxy);

CV=SqrtSumCovXY/meanSumXY;
CVbis=std(sumxy)/meanSumXY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZE DATA, (X) and (Y) INDEPENDENTLY, BEFORE STD CALCULATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx=x/mean(x);
yy=y/mean(y);
ssumxy=xx+yy;

vvarx=var(xx);
vvary=var(yy);
ccovxy=cov(xx,yy);
%............
ssumCovXY=sum(ccovxy(:));
SSqrtSumCovXY=sqrt(ssumCovXY); %Output.
%............

CCV=SSqrtSumCovXY;
CCVbis=std(ssumxy);

%....
SqrtSumCovXYout=[SqrtSumCovXY,SSqrtSumCovXY]
CVout=[CV,CVbis]
CCVout=[CCV,CCVbis]

%%%%%%
%PLOT:
%%%%%%
figure(1)
subplot(3,2,1)
plot(time,x1,time,x2,time,x3)
subplot(3,2,2)
plot(time,y1,time,y2,time,y3)
subplot(3,2,3)
plot(time,x)
subplot(3,2,4)
plot(time,y)
subplot(3,2,5)
plot(time,xx)
axis([0 +inf, 0 2])
subplot(3,2,6)
plot(time,yy)
axis([0 +inf, 0 2])

%%%%%%%%%%%%%%%%%%
%USING MY PROGRAM:
%%%%%%%%%%%%%%%%%%
%................................
[STDcase1,starSTDcase1]=mySTDforSeveralVariables([x1,x2,x3,y1,y2,y3]);
%................................
[STDcase2,starSTDcase2]=mySTDforSeveralVariables([x,y]);
%................................
[STDcase1phy,starSTDcase1phy]=mySTDforSeveralVariables([x1,x2,x3]);
[STDcase1zoo,starSTDcase1zoo]=mySTDforSeveralVariables([y1,y2,y3]);
%................................
[STDcase2phy,starSTDcase2phy]=mySTDforSeveralVariables(x);
[STDcase2zoo,starSTDcase2zoo]=mySTDforSeveralVariables(y);
%................................

disp('=================')
STDcase1case2=[STDcase1,STDcase2]
STDcase1case2phy=[STDcase1phy,STDcase2phy]
STDcase1case2zoo=[STDcase1zoo,STDcase2zoo]
disp('=================')
starSTDcase1case2=[starSTDcase1,starSTDcase2]
starSTDcase1case2phy=[starSTDcase1phy,starSTDcase2phy]
starSTDcase1case2zoo=[starSTDcase1zoo,starSTDcase2zoo]
disp('=================')

