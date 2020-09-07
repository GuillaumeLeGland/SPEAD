function [iSdin,iDrate,fsin,ttime,Wrad] = jamstecrest_funsin1D(nfreqs,ndays,t0,deltat,keyNutrientPulses,keyNutrientSupplyFrequencyConstant)

%===================================================================================
% $$$ %...................................................................................
% $$$ PulseFrequencyMax = 64;
% $$$ PulseFrequencyMin =  1;
% $$$ %...................................................................................
% $$$ % $$$ PulseFrequencyDel = (PulseFrequencyMax - PulseFrequencyMin) / (nfreqs - 1);
% $$$ % $$$ PulseFrequencyRng = [PulseFrequencyMin:PulseFrequencyDel:PulseFrequencyMax];
% $$$ %...................................................................................
% $$$ % $$$ PulseFrequencyDel = (log2(PulseFrequencyMax) - log2(PulseFrequencyMin)) / (nfreqs - 1);
% $$$ % $$$ PulseFrequencyRng = 2.^[log2(PulseFrequencyMin):PulseFrequencyDel:log2(PulseFrequencyMax)];
% $$$ %...................................................................................
% $$$ % $$$ %%PulseFrequency = [1,2,4,8,16,32,64];
% $$$ % $$$ PulseFrequency = [1,2,4,8,16,32,64,128];
% $$$ % $$$ PulseFrequency = [1,2,4,8,16,32,64,128,256,512];
% $$$ % $$$ PulseFrequency = 2.^([-3:8]);
% $$$ % $$$ PulseFrequency = 2.^([-2:9]);
% $$$ PulseFrequency = 2.^([-1:8]);
% $$$ %...................................................................................
% $$$ %%PulseFrequencyDel = (length(PulseFrequency) - 1) / (nfreqs - 1);
% $$$ %%PulseFrequencyRng = interp1([1:length(PulseFrequency)],PulseFrequency,[1:PulseFrequencyDel:length(PulseFrequency)])
% $$$ %...................................................................................
% $$$ PulseFrequencyRng = PulseFrequency; 
% $$$ %...................................................................................
% $$$ PulsePeriodStar = fliplr(PulseFrequencyRng/max(PulseFrequencyRng))';
% $$$ %...................................................................................
%===================================================================================
%...................................................................................
if     strcmp(keyNutrientSupplyFrequencyConstant,'yes');
    %...............................................................................
    %%PulseFrequency = 1*ones(nfreqs,1); 
    PulseFrequency = 4*ones(nfreqs,1); 
    %...............................................................................
elseif strcmp(keyNutrientSupplyFrequencyConstant,'not');
    %...............................................................................
    DoublingSeries = [1,2,4,8,16,32,64,128,256,512];
    %...............................................................................
    PulseFrequency = DoublingSeries(1:nfreqs)';
    %...............................................................................
end
%...................................................................................
PulsePeriodStar = 1./PulseFrequency; %From [0 - 1] 
%...................................................................................
%===================================================================================
%DEFINE THE PERIOD OF THE PERTURBATION: 
%...................................................................................
%%deltatpulse = 1;
deltatpulse = deltat;
%...................................................................................
ntyears = 3;
ttmax = ndays*ntyears; %I will cut out the first and third year afterwards.
time = [deltatpulse:deltatpulse:ttmax]; 
Time = ones(nfreqs,1)*time;
ntime = length(time); 
%...................................................................................
%%Tperiod = ndays; %(for SST seasonality)
%%Tperiod = ndays/4; %(for DIN supply)
%%Tperiod = ndays/30; %(for DIN supply)
%...................................................................................
Tperiod = ndays*PulsePeriodStar; 
%...................................................................................
%%tau = 0;
tau = Tperiod/4; 
Tau = tau*ones(1,ntime);
%...................................................................................
ttime = Time-Tau;
wrad = (2*pi)./Tperiod;
Wrad = wrad*ones(1,ntime);
%...................................................................................
%===================================================================================
%-----------------------------------------------------------------------------------
%NOTE: To make fsin to be between [0 - 1] use "Amp = 0.5" with "ysin = (Amp*sin(Wrad.*ttime) + Amp)".
%-----------------------------------------------------------------------------------
%...................................................................................
fsin = zeros(nfreqs,ntime);
%...................................................................................
% $$$ %%Amp = eps; %Amplitude [%] (for constant values)
% $$$ %%Amp = (1/3); %Amplitude [%] (for SST seasonality)
% $$$ Amp = 0.5; %Amplitude [%] (for DIN supply)
%%Amp = 0.6; %Amplitude [%] (for DIN supply) %original.
%%Amp = 0.9; %Amplitude [%] (for DIN supply)
Amp = 1.0; %Amplitude [%] (for DIN supply)
%...................................................................................
%%npower = 1.0;
npower = 2.0;
%...................................................................................
% $$$ ysin = (Amp*sin(Wrad.*ttime) + 1.0); %Sinusoidal variable centered at 1.0 [n.d.]
%...................................................................................
ysin = (0.5*sin(Wrad.*ttime) + 0.5); %Sinusoidal variable between [0 - 1] [n.d.] 
%...................................................................................
ymax = max(ysin,[],2)*ones(1,ntime);
ymin = min(ysin,[],2)*ones(1,ntime);
%...................................................................................
xsin = (ysin-ymin).^npower ./ (ymax-ymin).^npower; %Sinusoidal variable between [0 - 1] [n.d.] 
%...................................................................................
% $$$ fsin = ymin + 2.0*(Amp)*xsin.^npower; %Sinusoidal pulsed function centered at 1.0 [n.d.]
%...................................................................................
fsin = xsin.^npower; %Sinusoidal pulsed function between [0 - 1] [n.d.]
%...................................................................................
time1 = (ntime/ntyears)*1 + 1;
time2 = (ntime/ntyears)*2;
%...................................................................................
ymax = ymax(:,[time1:time2]); %Use only middle year.
ymin = ymin(:,[time1:time2]); %Use only middle year.
%...................................................................................
ysin = ysin(:,[time1:time2]); %Use only middle year.
xsin = xsin(:,[time1:time2]); %Use only middle year.
fsin = fsin(:,[time1:time2]); %Use only middle year.
%...................................................................................
% $$$ figure(1)
% $$$ subplot(2,2,1)
% $$$ imagesc(ysin)
% $$$ colorbar
% $$$ title('ysin')
% $$$ grid on
% $$$ subplot(2,2,2)
% $$$ imagesc(xsin)
% $$$ colorbar
% $$$ title('xsin')
% $$$ grid on
% $$$ subplot(2,2,4)
% $$$ imagesc(fsin)
% $$$ colorbar
% $$$ title('fsin')
% $$$ grid on
% $$$ pause(1)
% $$$ %%return
%...................................................................................
%===================================================================================
if strcmp(keyNutrientPulses,'Continuous')
%...................................................................................
fsincont = fsin; 
%...................................................................................
%===================================================================================
elseif strcmp(keyNutrientPulses,'Discrete')
%...................................................................................
ysin = (0.5*sin(Wrad.*ttime) + 0.5); %Sinusoidal variable between [0 - 1] [n.d.] 
%...................................................................................
for jdepth = 1:nfreqs
    %%Jpulses = [0:(Tperiod(jdepth)/deltat):(ttmax/deltat)];
    Jpulses = find(ysin(jdepth,:) == 1.0); %USAR ESTE MEJOR.
    Jpulses(1) = t0/deltat; %Put first pulse at t0
    Jpulses = Jpulses(2:end); %Remove first pulse at t0
    fsin(jdepth,Jpulses) = 1.0; 
end
%...................................................................................
time1 = (ntime/ntyears)*1 + 1;
time2 = (ntime/ntyears)*2;
%...................................................................................
ysin = ysin(:,[time1:time2]); %Use only middle year.
fsin = fsin(:,[time1:time2]); %Use only middle year.
%...................................................................................
intysin = sum(ysin*deltat,2);
intfsin = sum(fsin*deltat,2);
%...................................................................................
arearatio = intysin./intfsin; 
%...................................................................................
fsindisc = fsin; 
%...................................................................................
end %endif 
%===================================================================================
%...................................................................................
%HAVING DIN FLUX (CAN BE CONSTANT OR SEASONAL) AND OVERFLOW DILUTION:
%...................................................................................
Drate0  =  0.10; %Constant dilution specific rate [d-1]
dinext0 = 20.00; %Concentration of DIN exterior to my volume of water [mmolN*m-3]
%...................................................................................
minDrate0  =  0.1; %[d-1] 
minDINext0 =  1.0; %[mmolN*m-3] 
%...................................................................................
maxDrate0  =     minDrate0; %[d-1] 
maxDINext0 = 100*minDINext0; %[mmolN*m-3] 
%...................................................................................
%%DINloadperyear0 = 0.69*ndays; %[mmolN*m-3*d-1] x [days] = [mmolN*m-3] 
DINloadperyear0 = 0.54688*ndays; %[mmolN*m-3*d-1] x [days] = [mmolN*m-3] 
%...................................................................................
dinext0 = dinext0 / 2; %Decrease by a half to avoid negative concentrations during the ode45
DINloadperyear0 = DINloadperyear0 / 2;
%...................................................................................
DINfluxperyear0 = DINloadperyear0 + (minDINext0 * Drate0)*ndays; %[mmolN*m-3] 
%...................................................................................
%===================================================================================
if     strcmp(keyNutrientPulses,'Continuous')
    %...............................................................................
    DINext0 = dinext0*ones(size(fsin)); %make it same size than fsin.
    %...............................................................................
elseif strcmp(keyNutrientPulses,'Discrete')
    %...............................................................................
    for jdepth = 1:nfreqs
	npulses = length(find(fsin(jdepth,:) == 1.0));
	dinpulse0 = DINloadperyear0 / npulses; %Concentration added per discrete pulse [mmolN*m-3]
	%%dinpulse0 = DINloadperyear0 * PulsePeriodStar(jdepth); %Concentration added per discrete pulse [mmolN*m-3]
	zdinext0(jdepth,:) = dinpulse0 / (Drate0 * deltat); %External concentration maximum of DIN [mmolN*m-3]
    end
    %...............................................................................
    DINext0 = zdinext0 * ones(1,size(fsin,2)); %make it same size than fsin.
    %...............................................................................
end %endif 
%........................................................................
%===================================================================================
%PULSED DILUTION RATE (WITH CONSTANT DIN EXTERIOR):
%...................................................................................
% $$$ % $$$ iDrate =  Drate0  * fsin; %Sinusoidal dilution specific rate [d-1]
% $$$ % $$$ iDINext = DINext0 * ones(size(fsin));
% $$$ %...................................................................................
% $$$ % $$$ iDrate =  maxDrate0  * fsin; %Sinusoidal dilution specific rate [d-1]
% $$$ % $$$ iDINext = maxDINext0 * ones(size(fsin));
% $$$ %...................................................................................
% $$$ iDrate =   minDrate0 + (maxDrate0  * fsin); %Sinusoidal dilution specific rate [d-1]
% $$$ iDINext = minDINext0 + (maxDINext0 * ones(size(fsin)));
%...................................................................................
%===================================================================================
%PULSED DIN EXTERIOR CONCENTRATION (WITH CONSTANT DILUTION):
%...................................................................................
% $$$ iDrate  = Drate0  * ones(size(fsin));
% $$$ iDINext = DINext0 .* fsin; %Sinusoidal concentration of DIN exterior to my volume of water [mmolN*m-3]
%...................................................................................
% $$$ %%iDrate  = maxDrate0  * ones(size(fsin));
% $$$ %%iDINext = maxDINext0 .* fsin; %Sinusoidal concentration of DIN exterior to my volume of water [mmolN*m-3]
%...................................................................................
% $$$ iDrate  = minDrate0  + (maxDrate0  * ones(size(fsin)));
% $$$ iDINext = minDINext0 + (maxDINext0 .* fsin); %Sinusoidal concentration of DIN exterior to my volume of water [mmolN*m-3]
%...................................................................................
iDrate  = Drate0 * ones(size(fsin)); %Constant dilution rate [d-1] 
iDINext = minDINext0 + (DINext0 .* fsin); %Sinusoidal concentration of DIN exterior to my volume of water [mmolN*m-3]
%...................................................................................
%===================================================================================
%PULSED BOTH DILUTION RATE AND DIN EXTERIOR CONCENTRATION: 
%...................................................................................
% $$$ % $$$ iDrate =  Drate0 + (Drate0 * fsin); %Sinusoidal dilution specific rate [d-1]
% $$$ % $$$ iDINext = DINext0 .* fsin; %Sinusoidal concentration of DIN exterior to my volume of water [mmolN*m-3]
% $$$ %...................................................................................
% $$$ % $$$ iDrate =  maxDrate0  * fsin; %Sinusoidal dilution specific rate [d-1]
% $$$ % $$$ iDINext = maxDINext0 .* fsin; %Sinusoidal concentration of DIN exterior to my volume of water [mmolN*m-3]
% $$$ %...................................................................................
% $$$ iDrate =  minDrate0  + (maxDrate0  * fsin); %Sinusoidal dilution specific rate [d-1]
% $$$ iDINext = minDINext0 + (maxDINext0 .* fsin); %Sinusoidal concentration of DIN exterior to my volume of water [mmolN*m-3]
%...................................................................................
%===================================================================================
%CONSTANT BOTH DILUTION AND DIN EXTERIOR CONCENTRATION: 
%...................................................................................
% $$$ iDrate =   Drate0 * ones(size(fsin)); %Constant dilution specific rate [d-1]
% $$$ iDINext = DINext0 * ones(size(fsin)); %Constant concentration of DIN exterior to my volume of water [mmolN*m-3]
% $$$ %...................................................................................
% $$$ % $$$ iDrate =   minDrate0 * ones(size(fsin)); %Constant dilution specific rate [d-1]
% $$$ % $$$ iDINext = minDINext0 * ones(size(fsin)); %Constant concentration of DIN exterior to my volume of water [mmolN*m-3]
% $$$ %...................................................................................
% $$$ % $$$ iDrate =   maxDrate0 * ones(size(fsin)); %Constant dilution specific rate [d-1]
% $$$ % $$$ iDINext = maxDINext0 * ones(size(fsin)); %Constant concentration of DIN exterior to my volume of water [mmolN*m-3]
% $$$ %...................................................................................
% $$$ % $$$ iDrate =  minDrate0  + (maxDrate0  * ones(size(fsin))); %Constant dilution specific rate [d-1]
% $$$ % $$$ iDINext = minDINext0 + (maxDINext0 * ones(size(fsin))); %Constant concentration of DIN exterior to my volume of water [mmolN*m-3]
%...................................................................................
%===================================================================================
%...................................................................................
iSdin =  iDrate .* iDINext; %Flux of nutrient supply [mmolN*m-3*d-1]
%...................................................................................
%===================================================================================
%...................................................................................
iMdin = zeros(size(fsin)); %Zero nutrient losses here (I will add them later as fraction of plankton, DIN and PON concentrations).
%...................................................................................
intFsin = sum(fsin*deltatpulse,2); 
%...................................................................................
intSdin = sum(iSdin*deltatpulse,2); %Integrated flux of nutrient supply for the whole run [mmolN*m-3]
%...................................................................................
% $$$ iMdin = (intSdin/ttmax)*ones(size(fsin)); %Constant nutrient losses [mmolN*m-3*d-1]
%...................................................................................
intMdin = sum(iMdin*deltatpulse,2);
%...................................................................................
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%...................................................................................
%test: Check that total nutrient supply is constant for different frequencies: 
%...................................................................................
diffintSdin = diff(intSdin); 
Inonzero = find(abs(diffintSdin) > sqrt(eps));
if length(Inonzero) > 0
    intSdin 
    disp('Warning!!! The total nutrient supply should be *constant* for different frequencies!!!')
    pause 
end
%...................................................................................
intSdin
DINfluxperyear0
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%===================================================================================
% $$$ %...................................................................................
% $$$ %Check: 
% $$$ if abs(intSdin - intMdin) < 1d-6 %Almost zero. 
% $$$     display('Okay')
% $$$ else
% $$$     display('Error!!! DIN Supply and Losses *must* be balanced on a yearly integrated basis')
% $$$ end
%...................................................................................
figure(10)
subplot(2,2,1)
imagesc(iDINext,[0-eps max(iDINext(:))]) 
colorbar 
title('DIN external [mmolN*m-3]')
grid on
subplot(2,2,2)
imagesc(iDrate,[0-eps max(iDrate(:))]) 
colorbar 
title('Dilution rate [d-1]')
grid on
subplot(2,2,3) 
imagesc(iSdin,[0-eps max(iSdin(:))])
colorbar 
title('Supply of DIN [mmolN*m-3*d-1]')
grid on
subplot(2,2,4) 
imagesc(iMdin,[0-eps max(iMdin(:))]) 
colorbar 
title('Losses of DIN [mmolN*m-3*d-1]')
grid on
%...................................................................................
print('-dpng','-r300','jamstecrest_gaussecomodel1D_fig005.png')
%...................................................................................
%%return 
pause(1)
close all 
%...................................................................................
%===================================================================================
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%...................................................................................
% $$$ Mdin = zeros(1,ntime/ntyears); %Remove nutrient losses (just for checking things)
%...................................................................................
% $$$ Sdin = zeros(1,ntime/ntyears); %Remove nutrient supply (just for checking things)
%...................................................................................
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
return

