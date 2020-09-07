function [aveSDIN,aveNTOT,aveMUP,aveMUZ,aveFPHYT,aveFZOO,avePHYT,aveZOO,aveDIN,avePON,avelogESDphyAve,avelogESDphyStd,aveESDphyAve,aveESDphyCV,aveSindex,aveEindex,aveShannonEntropyTheoretical] = jamstecrest_loopruns_longtermaverage(SDIN,NTOT,MUP,MUZ,FPHYT,FZOO,PHYT,ZOO,DIN,PON,logESDphyAve,logESDphyStd,ESDphyAve,ESDphyCV,Sindex,Eindex,ShannonEntropyTheoretical)

%====================================================================
%....................................................................
aveSDIN = squeeze(mean(SDIN,2)); %[depth,frequency]
aveNTOT = squeeze(mean(NTOT,2)); %[depth,frequency]
%....................................................................
aveMUP = squeeze(mean(MUP,2)); %[depth,frequency] 
aveMUZ = squeeze(mean(MUZ,2)); %[depth,frequency]
%....................................................................
aveFPHYT = squeeze(mean(FPHYT,2)); %[depth,frequency]
aveFZOO = squeeze(mean(FZOO,2)); %[depth,frequency]
%....................................................................
avePHYT = squeeze(mean(PHYT,2)); %[depth,frequency]
aveZOO = squeeze(mean(ZOO,2)); %[depth,frequency]
aveDIN = squeeze(mean(DIN,2)); %[depth,frequency]
avePON = squeeze(mean(PON,2)); %[depth,frequency]
%....................................................................
avelogESDphyAve = squeeze(mean(logESDphyAve,2)); %[depth,frequency]
avelogESDphyStd = squeeze(mean(logESDphyStd,2)); %[depth,frequency]
aveESDphyAve = squeeze(mean(ESDphyAve,2)); %[depth,frequency]
aveESDphyCV  = squeeze(mean(ESDphyCV,2)); %[depth,frequency]
%....................................................................
aveShannonEntropyTheoretical = squeeze(mean(ShannonEntropyTheoretical,2)); %[depth,frequency]
aveSindex = squeeze(mean(Sindex,2)); %[depth,frequency]
aveEindex = squeeze(mean(Eindex,2)); %[depth,frequency]
%....................................................................
%====================================================================
%********************************************************************
return





