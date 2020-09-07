function [muKMCchl,muKMCchlTheoretical] = myShannonContinousEntrophyFromLanSmith(simchl,logESD,logESDave,logESDstd);

%====================================================================
%....................................................................
simchl = simchl(:);
%....................................................................
logESD = logESD(:);
%....................................................................
%====================================================================
%calculate the discrete Shannon Index based on chl: 
%....................................................................
minChl = 0; %How much should I use? 
TotChl = sum(simchl); 
Npft = length(simchl); 
%....................................................................
HindexChl = 0; 
for j = 1:Npft
    if(simchl(j) > minChl ) 
	pESDchlj = simchl(j) / TotChl; 
	HindexChl = HindexChl - pESDchlj * log(pESDchlj); 
    end
end
%....................................................................
%====================================================================
%Calculate the Shannon Entropy, mu [see Quintana et al. 2008] 
%expected value of Y = logESD (here denoted EY) 
%and expected value of Y*Y (here denoted EY2). 
%....................................................................
mulESDchl = sum(simchl .* logESD) / TotChl; 
Ynormchl = logESD - mulESDchl; 
%....................................................................
%====================================================================
%....................................................................
EYchl = 0.0d0;
EYchl2 = 0.0d0;
for j = 1:Npft
    pESDchlj =  simchl(j)   / TotChl;
    EYchl   = EYchl  + pESDchlj * Ynormchl(j);
    EYchl2  = EYchl2 + pESDchlj * Ynormchl(j)^2;
end
%....................................................................
%====================================================================
%(Quintana et al. 2008, p. 77, below eqn. 6)
%the minimum value (0.20) was chosen by trial and error to avoid NaNs for muKMC 
%....................................................................
sKernchl = 1.06 * sqrt(EYchl2 - EYchl^2) / Npft^0.20d0; 
%....................................................................
%%sKernchl = max( [2.0d-1 , 1.06 * sqrt(EYchl2 - EYchl^2) / Npft^0.20d0] ); 
%....................................................................
%====================================================================
% Continuous Shannon Entropy (eq. 9 of Quintana et al. 2008). 
% Here using concentrations instead of counts;  i. e., 
% summing over the uniquely defined PFTs, weighted by the concentrations, 
% instead of summing over the number of samples counted.  
% This is mathematically consistent.
%....................................................................
muKMCchl = EYchl + log(sKernchl*sqrt(2*pi));
%....................................................................
for k = 1:Npft
    pESDchlk =    simchl(k) / TotChl; 
    sum1chl = 0.0d0; 
    for j = 1:Npft
	pESDchlj = simchl(j) / TotChl;
	sum1chl = sum1chl + pESDchlj*exp( -( ((Ynormchl(k)-Ynormchl(j))/sKernchl)^2 )/2 );
    end
    muKMCchl = muKMCchl - pESDchlk*log(sum1chl);
end
%....................................................................
%====================================================================
%EXACT SOLUTION FOR THE CONTINOUS LOG-NORMAL DISTRIBUTION: 
%....................................................................
logESDstdUnbias = logESDstd*sqrt(Npft/(Npft-1)); 
%....................................................................
muKMCchlTheoretical001 = (1/2) + log(logESDstdUnbias * sqrt(2*pi)) + logESDave; 
%....................................................................
%====================================================================
%FROM WIKIPEDIA: 
% <https://en.wikipedia.org/wiki/Normal_distribution> 
% <http://www.biopsychology.org/norwich/isp/chap8.pdf> 
%....................................................................
muKMCchlTheoretical002 = (1/2) * log((2*pi)*(logESDstd.^2)) + (1/2);
%....................................................................
%====================================================================
%....................................................................
%%muKMCchlTheoretical = muKMCchlTheoretical001; %WEIRD...
muKMCchlTheoretical = muKMCchlTheoretical002; %MAYBE OKAY.
%....................................................................
%********************************************************************
return

