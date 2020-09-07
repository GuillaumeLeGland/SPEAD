function [RandomBinOut]=myRandomBin(DistributionType,nptos,upperlimit,deltabin)

%***********************************************************************
%Use: [RandomBinOut]=myRandomBin(DistributionType,nptos,upperlimit,deltabin)
%
%If RandomBin > 1 the paramters increases a percentage: abs(RandomBin-1)
%If RandomBin < 1 the paramters decreases a percentage: abs(RandomBin-1)
%***********************************************************************
%.......................................................................
if strcmp(DistributionType,'Uniforme')
    PDFname = 'Uniform';
elseif strcmp(DistributionType,'Gaussian')
    PDFname = 'Norm';
end
%.......................................................................
if deltabin > 0
    %========================================
    %GENERATE RANDOM NUMBERS BETWEEN 0 AND 1:
    %========================================
    if strcmp(DistributionType,'Uniform')
	%............
	RandomDistrib = rand(1,nptos)*upperlimit; %Random value between Zero and Upperlimit.
	%............
    elseif strcmp(DistributionType,'Gaussian')
	%............
	Param1=0.5; %mean.
	Param2=0.2; %sigma.
	msize=1;
	nsize=nptos;
	%............
	RandomDistrib = random(PDFname,Param1,Param2,msize,nsize);
	%............
	I=find(RandomDistrib>1);
	J=find(RandomDistrib<0);
	RandomDistrib(I)=0.9;
	RandomDistrib(J)=0.1;
	%............
    end
    %.........................

    %=====================
    %OBTAIN DISCRETE BINS:
    %=====================
    RandDiscrete=ceil(RandomDistrib/deltabin)*deltabin; %Bin the values in ranges of 20% [0-0.2, 0.2-0.4, 0.4-0.6, 0.6-0.8, 0.8-1.0]
    RandDiscreteMiddle=RandDiscrete-(0.5*deltabin); %Use the middel point of each range %[0.1, 0.3, 0.5, 0.7, 0.9]
    RandDiscreteMiddleCenteredAtOne=(1-upperlimit/2)+RandDiscreteMiddle; %[0.6, 0.8, 1.0, 1.2, 1.4]

    %======
    %OUPUT:
    %======
    RandomBin=RandDiscreteMiddleCenteredAtOne; %[0.6, 0.8, 1.0, 1.2, 1.4]
else
    RandomBin=ones(1,nptos);
end
%..........................

%OUTPUT:
%.................................
RandomBinOut(1:nptos)=RandomBin; %USAR ESTE.
%.................................
% $$$ SmallNoise = (rand(1)-0.5)*sqrt(eps);
% $$$ RandomBinOut(1:nptos) = RandomBin + SmallNoise;
%.................................

return
%***************************************
%..................
nptos=1000;
pcntmax=0.5;
%..................
%%binsize=0.05;
binsize=0.1;
%%binsize=0.2;
%..................
upperlimit=pcntmax*2;
%..................
deltabin=binsize*upperlimit; %For values between 0 and Upperlimit.
%..................
nparamsX=2;
% $$$ nparamsX=6;
%..................
%%DistributionType='Uniform';
DistributionType='Normal';
%..................
RandomBinX=[];
for i=1:nptos
    [RandomBinXi]=modelNPHC_RandomBin2(DistributionType,nparamsX,upperlimit,deltabin);
    RandomBinX=[RandomBinX;RandomBinXi];
end

xRandomBinX=RandomBinX;
I=find(RandomBinX==1);
xRandomBinX(I)=nan;
figure(1)
hist(xRandomBinX(:),100)
axis([0 2, 0 nptos/2])
