%function [jday,jjday,newday] = jamstecrest_daycounter(jday,iTime)
% jjday is a useless variable, always equal to jday
function [jday,newday] = jamstecrest_daycounter(jday,iTime) 

% jday0 = jday; %Previous day.
% jday = ceil(iTime); %Current day.
% jday1 = jday; %Current day.
% if jday0 == jday1
%     newday = 'not';
%     jjday  = jday0; %Previous day.
% else
%     newday = 'yes';
%     jjday  = jday1; %Current day.
% end

jdaynew = ceil(iTime); % Current day
if jdaynew == jday
    newday = 'not';
else
    newday = 'yes';
    jday = jdaynew;
end

return
    
