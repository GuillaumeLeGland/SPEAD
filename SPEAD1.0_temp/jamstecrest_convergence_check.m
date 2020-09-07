function [conv_P,conv_xave,conv_yave,conv_xxvar,conv_yyvar,conv_xycor] = jamstecrest_convergence_check(P,xave,yave,xxvar,yyvar,xycov,ndays)
% Check if the simulation is convergent (Le Gland, 17/07/2020) using
% variations between years

% Proportional variation compared with the previous year (total phytoplankton)
conv_P     = (P(:,ndays+1:end)-P(:,1:end-ndays)) ./ P(:,1:end-ndays);
% Variation compared with the previous year normalized by the range of trait values (xave and yave)
conv_xave  = (xave(:,ndays+1:end)-xave(:,1:end-ndays)) ./ (max(xave(:))-min(xave(:)));
conv_yave  = (yave(:,ndays+1:end)-yave(:,1:end-ndays)) ./ (max(yave(:))-min(yave(:)));
% Proportional variation compared with the previous year (variances)
conv_xxvar = (xxvar(:,ndays+1:end)-xxvar(:,1:end-ndays)) ./ xxvar(:,1:end-ndays);
conv_yyvar = (yyvar(:,ndays+1:end)-yyvar(:,1:end-ndays)) ./ yyvar(:,1:end-ndays);
% variation of correlation compared with previous year (no unit)
xycor = xycov./sqrt(xxvar.*yyvar);
conv_xycor = xycor(:,ndays+1:end)-xycor(:,1:end-ndays);

disp('largest difference between year n and n-1 on P');
max(max(abs(conv_P(:,end-ndays+1:end))))
disp('largest difference between year n and n-1 on xave');
max(max(abs(conv_xave(:,end-ndays+1:end))))
disp('largest difference between year n and n-1 on yave');
max(max(abs(conv_yave(:,end-ndays+1:end))))
disp('largest difference between year n and n-1 on xxvar');
max(max(abs(conv_xxvar(:,end-ndays+1:end))))
disp('largest difference between year n and n-1 on yyvar');
max(max(abs(conv_yyvar(:,end-ndays+1:end))))
disp('largest difference between year n and n-1 on xycor');
max(max(abs(conv_xycor(:,end-ndays+1:end))))


end

