function [SqrtSumVarX,SqrtSumCovX,nSqrtSumVarX,nSqrtSumCovX,nvars]=mySTDforSeveralVariables(type,varargin);

%**********************************
%PROGRAMA "mySTDforSeveralVariables.m":
%How to calculate the Standard Deviation of the sum of several variables.
%
%Use:
%    
%1st way: [SqrtSumCovX,SqrtSumVarX] = mySTDforSeveralVariables(X)
%
%where each X has got to be an array where each row is an observation
%(eg. a time-series) and each column is a variable.
%
%2nd way: [SqrtSumCovX,SqrtSumVarX] = mySTDforSeveralVariables(x1(:),x2(:),...,xn(:))
%
%where each xi has got to be a column vector (eg. a time-series) of same size.
%
%
%Input:
%
% X (time,vari) array.
%
%Outputs:
%
%a) SqrtSumVarX:  Sqrt of the sum of Variances of X                     %(no del todo correcto, pues no incluye covarianzas)
%b) SqrtSumCovX:  Sqrt of the sum of Variances of X + Covariances of X  %(correcto)
%
%c) nSqrtSumVarX: Sqrt of the sum of Variances of X                    / number of time-series(*)
%d) nSqrtSumCovX: Sqrt of the sum of Variances of X + Covariances of X / number of time-series(*)
%
%(see also "egSTDforSeveralVariables.m" and 
%"THOMPSON_2008_#The_standard_deviation_of_the_sum_of_several_variables#_AMC_30.pdf").
%"TILMAN_1999_#The_ecological_consequences_of_changes_in_biodiversity__A_search_for_general principles#_Ecology_80.pdf
%
%(*) NOTA: Dividir por el numero de time-series "n" permite standarizar
%el resultado de manera que no dependa del numero de time-series
%utilizado (si no, a mayor numero de time-series, mayor valor de la STD).
%De esta forma puedo compara la STD obtenida con una sola curva temporal
%vs. la STD obtenida con muchas curvas temporales!!.
%**********************************

nvarargin=length(varargin);

if nvarargin==1
    Xinput=varargin{1};
elseif nvarargin>=2
    for i=1:nvarargin
	Xinput(:,i)=varargin{i};
    end
end

%................................
%--------------------------------
%NOTE:
%Use "DataRaw" in the case the input variable are Anomalies (ie. dimensionless values with mean equal to ~zero).
%Use "DataNormalizeByMeanPopulations" in the rest of the cases (ie. values with dimension and non-zero mean).
%--------------------------------
[mtimes,nvars] = size(Xinput); %(time,vars) %Make sure that X comes in this way!!!
%................................
if nvars == 0
    SqrtSumVarX  = nan;
    SqrtSumCovX  = nan;
    nSqrtSumVarX = nan;
    nSqrtSumCovX = nan;
    return
end
%................................
if strcmp(type,'DataRaw')
    %=========
    %Raw data:
    %=========
    Xnorm = Xinput;
elseif strcmp(type,'DataNormalizedByMeanTotal')
    %======================================================
    %Normalize data by total mean to make them dimensionless:
    %======================================================
    XT = mynansum(Xinput,2); %size(time,1)
    XTave = nanmean(XT);
    %%pause 
    for jvar=1:nvars
	%...............
	Xj = Xinput(:,jvar);
	%...............
	Xnorm(:,jvar) = Xj/XTave;
	%...............
    end
elseif strcmp(type,'DataNormalizedByMeanPopulations')
    %======================================================
    %Normalize data by each time-seris mean to make them dimensionless:
    %======================================================
    for jvar=1:nvars
	%...............
	Xj = Xinput(:,jvar);
	%...............
	meanXj = mean(Xj);
	%...............
	Xnorm(:,jvar) = Xj/meanXj;
	%...............
    end
end
X = Xnorm;
%................................

%%%%%%%%%%%%%%%%%%%%%%%
%CHANGE NaNs FOR Zeros:
%%%%%%%%%%%%%%%%%%%%%%%
Inan = find(isnan(X)==1);
X(Inan) = 0;

%================================================================
%a) Taking into account the covariances (in theory more correct): 
%sum of the variances plus sum of the covariances.
%================================================================
%--------------------------------------------------------------------
%NOTE: When the variables are strongly anti-correlated using this one
%does not always work (I can get 'zero' as a result).
%--------------------------------------------------------------------
covX = cov(X);
sumCovX = sum(covX(:)); %var(Phy + Zoo) = var(Phy) + var(Zoo) + 2*cov(Phy,Zoo)
SqrtSumCovX = sqrt(sumCovX); %Output.
SqrtSumCovX = real(SqrtSumCovX); %in case SqrtSumCovX is a negative value...

%========================================================================
%test: Old-school way; First sum up the variables and then calculate std:
%========================================================================
sumX = mynansum(X,2); %size(time,1)
SqrtSumCovXbis = std(sumX);
xdist = abs(SqrtSumCovX-SqrtSumCovXbis);

if xdist > 1d-10 | isfinite(xdist) == 0 %If "xdist" is too big or NaN.
    [SqrtSumCovX,SqrtSumCovXbis]; %should be the same.
    %............
    covX
    sumCovX
    %............
    SqrtSumCovX
    SqrtSumCovXbis
    %............
    xdist
    %............
    disp('SqrtSumCovX and SqrtSumCovXbis should be the same!!')
    pause
end
%.............

%===================================================================
%b) Without taking into account the covariances (kind of incorrect): 
%sum of the variances alone.
%===================================================================
traceCovX = trace(covX); %equivalent to sum(var(X)).
%.....
stdX = std(X);
varX = var(X);
sumStdX = mynansum(stdX);
sumVarX = mynansum(varX);
if abs(traceCovX-sumVarX) > 1d-6
    disp('Error! traceCovX and sumVarX should be the same!')
    [sumVarX,traceCovX] %should be the same.
    pause 
end
%.....
SqrtTraceCovX = sqrt(traceCovX);
SqrtSumVarX = sqrt(sumVarX);

%????????????????????????
% $$$ SqrtTraceCovX=traceCovX;
% $$$ SqrtSumVarX=sumVarX;
%????????????????????????

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZE BY DIVIDING BY THE NUMBER OF TIME-SERIES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...........................
% $$$ nSqrtSumVarX = sqrt(sumVarX); %sum of variances.
% $$$ nSqrtSumCovX = sqrt(sumCovX); %sum of variances + sum of covariances.
%...........................
nSqrtSumVarX = sqrt(sumVarX/nvars); %sum of variances. %USAR ESTOS!!!!
nSqrtSumCovX = sqrt(sumCovX/nvars); %sum of variances + sum of covariances.
%...........................
% $$$ nSqrtSumVarX = SqrtSumVarX/nvars; %ESTA MAL!!!! (see example below with "varx")
% $$$ nSqrtSumCovX = SqrtSumCovX/nvars; 
%...........................

%...........................
% $$$ %test
% $$$ whos X 
% $$$ 
% $$$ stdX
% $$$ varX
% $$$ 
% $$$ sumStdX
% $$$ sumVarX
% $$$ 
% $$$ SqrtSumVarX 
% $$$ SqrtSumVarXok = sum(stdX) %creo que deben salir iguales!
% $$$ 
% $$$ varX
% $$$ sumVarX
% $$$ sqrt(sumVarX)
%...........................

return
%**************************************
%Copy-paste this bit in the command window to check the fiability of "mySTDforSeveralVariables.m"
%...........
x1 = rand(1,365);
x2 = x1;
x3 = x1;
x4 = x1;
%...........
% $$$ X1 = x1(:);
% $$$ X2 = [x1(:),x2(:),x3(:),x4(:)];
%...........
X2 = [x1(:),x2(:),x3(:),x4(:)]/4;
X1 = sum(X2,2);
%...........
figure(1)
subplot(2,2,1)
plot([1:365],X1)
grid on
subplot(2,2,2)
plot([1:365],X2)
grid on
%...........
%%type = 'DataRaw';
type = 'DataNormalizedByMeanPopulations';
type = 'DataNormalizedByMeanTotal';
%...........
[SqrtSumVarX1,SqrtSumCovX1,nSqrtSumVarX1,nSqrtSumCovX1] = mySTDforSeveralVariables(type,X1);
[SqrtSumVarX2,SqrtSumCovX2,nSqrtSumVarX2,nSqrtSumCovX2] = mySTDforSeveralVariables(type,X2);
%...........
SqrtSumVarX1
SqrtSumVarX2
%...........
nSqrtSumVarX1 %Han de salir iguales (for type = 'DataNormalizedByMeanPopulations')
nSqrtSumVarX2
%...........

%OJO, HAY ALGO MAL!!!
%Para 'type = 'DataRaw' deberia tener SqrtSumVarX1 = SqrtSumVarX2, pero
%obtengo que el primero es el doble que el segundo! Creo que el problema
%va mas o menos por lo siguiente:

varx = [0.1,0.1,0.1,0.1]

std1 = sqrt(sum(varx))
std2 = sum(sqrt(varx))

%No es lo mismo hacer "std1" que hacer "std2", y obtengo que uno es el doble del otro!!

