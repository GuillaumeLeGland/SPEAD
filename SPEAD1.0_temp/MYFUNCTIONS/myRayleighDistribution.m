function [FDF,PDF,PDFchap,PDFstats,PDFparms,xbins]=myRayleighDistribution(xdata,varargin)

nvarargin=length(varargin);    
if nvarargin==0 %Calculo simga from the data.
    %...................
    N=length(xdata);
    lambda = sqrt( (1/(2*N) * sum(xdata.^2)) );
    %...................
    [lambdaBis] = mytbxRaylfit(xdata);
    kpower=2; 
    %...................
    %test:
    xlambda=[lambda,lambdaBis] %should be the same!
    %...................
    %Weibull PDF:
    Izero=find(xdata==0);
    xxdata=xdata;
    xxdata(Izero)=eps;
    lambdaPDFrayl=mytbxRaylfit(xxdata);
    [w]=mytbxWblfit(xxdata);
    wlambda=w(1); %scale param.
    wkpower=w(2); %shape param.
    %...................
elseif nvarargin==1 %Impose lambda y kpower a-priori.
    lambda=varargin{1}; 
    kpower=varargin{2};
end

%?????????????????
% $$$ wkpower=1.5; %For NPHCpaper.
%?????????????????
%..................
wkpowerMin=1.25; %For NPHCpaper.
if wkpower < wkpowerMin
    wkpower =  wkpowerMin;
end
%..................

%%%%%%%%%%%%%%%%%%
%OBTENGO FDF REAL:
%%%%%%%%%%%%%%%%%%
%........................
nbins=100;
xmin=min(xdata(:));
xmax=max(xdata(:));
deltax=(xmax-xmin)/(nbins-1);
xbins=[xmin:deltax:xmax];
%........................
[FDF] = histc(xdata,xbins); %Frecuency Distrib Function.
%........................
dAreasFDF = deltax.*FDF;
AreaFDF = sum(dAreasFDF)
%........................

%====================
%STATISTICAL MOMENTS:
%====================
FDFmean = mean(xdata);
FDFvar = var(xdata); 
%....
FDFmean_FDFvar=[FDFmean,FDFvar]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO PDF (ie. AREA BAJO LA CURVA HA DE SER IGUAL A 1.0):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%........................
PDF = FDF/AreaFDF; %Probability Distrib Function.
%........................
dAreasPDF = deltax.*PDF;
AreaPDF = sum(dAreasPDF)
%........................

%%%%%%%%%%%%%%%%%%%%%
%OBTENGO PDF TEORICA:
%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------
%NOTA: Usando esto deberia obtener que Weibull y Rayleight 
%son PDF equivalentes; pero no me sale!!!
%a) wlambda=lambda; %Rayleigh sigma.
%b) wkpower=2;
%----------------------------------------------------------
xchap = xbins;

%===============
%EQUATION TERMS:
%===============
%....................
%Rayleigh:
alfa1 = (xchap/lambda^2);
beta1 = (-0.5*(xchap/lambda).^2);
%....................
%Weibull:
alfa2 = (wkpower/wlambda)*(xchap/wlambda).^(wkpower-1);
beta2 = -(xchap./wlambda).^wkpower;
%....................

%===============
%FINAL EQUATION:
%===============
% $$$ alfa=alfa1; %Rayleigh.
% $$$ beta=beta1;
%....................
alfa=alfa2; %Weibull.
beta=beta2;
%....................
PDFchap = alfa.*exp(beta);
%....................
dAreasPDFchap = deltax.*PDFchap;
AreaPDFchap = sum(dAreasPDFchap)
%........................

%====================
%STATISTICAL MOMENTS:
%====================
%........................
raylPDFmean = lambda*sqrt(pi/2);
raylPDFmedian = lambda*sqrt(log(4));
raylPDFmode = lambda;
raylPDFvar = 0.5*(4-pi)*lambda^2;
raylPDFstd = sqrt(raylPDFvar);
%........................
weibPDFmedian = wlambda*(log(2))^(1/wkpower);
weibPDFmode = wlambda*((wkpower-1)/wkpower)^(1/wkpower);
[weibPDFmean,weibPDFvar] = mytbxWblstat(wlambda,wkpower);
weibPDFstd = sqrt(weibPDFvar);
%........................
%test: 
m = wlambda .* gamma (1 + 1 ./ wkpower);
v = (wlambda .^ 2) .* gamma (1 + 2 ./ wkpower) - m .^ 2;
weibPDFmeanBis = weibPDFmean;
weibPDFvarBis = weibPDFvar;
if weibPDFmean ~= weibPDFmeanBis | weibPDFvar ~= weibPDFvarBis
    weibPDFmean,weibPDFmeanBis
    weibPDFvar,weibPDFvarBis
    disp('Error! They must have the same value!')
    pause
end
% $$$ [weibPDFmean,weibPDFvar]
% $$$ [weibPDFmeanBis,weibPDFvarBis]
% $$$ pause
% $$$ weibPDFmean = weibPDFmeanBis;
% $$$ weibPDFvar = weibPDFvarBis;
%........................

%........................
%Rayleigh:
% $$$ PDFmean = raylPDFmean;
% $$$ PDFmedian = raylPDFmedian;
% $$$ PDFmode = raylPDFmode;
% $$$ PDFvar = raylPDFvar;
% $$$ PDFstd = raylPDFstd;
% $$$ PDFlambda = lambda;
% $$$ PDFkpower = kpower;
%........................
%Weibull:
PDFmean = weibPDFmean;
PDFmedian = weibPDFmedian;
PDFmode = weibPDFmode;
PDFvar = weibPDFvar;
PDFstd = weibPDFstd;
PDFlambda = wlambda;
PDFkpower = wkpower;
%........................

%%%%%%%%%
%OUTPUTS:
%%%%%%%%%
%......................
PDFstats(1)=PDFmean;
PDFstats(2)=PDFmedian;
PDFstats(3)=PDFmode;
PDFstats(4)=PDFvar;
PDFstats(5)=PDFstd;
%......................
PDFparms(1)=PDFlambda;
PDFparms(2)=PDFkpower;
%......................
%....
PDFstats
PDFparms

%....................
% $$$ figure(10)
%....................
% $$$ subplot(2,2,1)
% $$$ hist(xdata,nbins);
% $$$ hold on
% $$$ plot(xbins,FDF,'r*')
% $$$ hold off
% $$$ %....................
% $$$ subplot(2,2,2)
% $$$ hp1=plot(xbins,PDF,'r*');
% $$$ hold on
% $$$ hp2=plot(xbins,PDFchap,'-b.')
% $$$ hold off
% $$$ set(hp1,'Color',[0.7 0.7 0.7])
% $$$ set(hp2,'Color',[0.1 0.0 0.0])
% $$$ set(gca,'Xlim',[0 1.5],'Ylim',[0 10.0])
% $$$ %....................
% $$$ subplot(2,2,3)
% $$$ plot(PDF,PDFchap,'*')
% $$$ xlabel('PDF')
% $$$ ylabel('PDFchap')
% $$$ %....................
% $$$ subplot(2,2,4)
% $$$ hb=bar(xbins,PDF,1.0,'grouped');
% $$$ hold on
% $$$ plot(xbins,PDFchap,'-r.')
% $$$ hold off
% $$$ xlabel('PDF')
% $$$ ylabel('PDFchap')
% $$$ set(gca,'Xlim',[0 1.5],'Ylim',[0 10.0])
% $$$ set(hb,'Facecolor',[0.7 0.7 0.7])
%....................
% $$$ pause
% $$$ close all
%....................

%**********************
return
figure(40)
plot(xbins,PDF,'-r+')
f = ezfit('gauss')

