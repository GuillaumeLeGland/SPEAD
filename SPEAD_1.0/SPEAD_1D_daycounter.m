function [jday,newday] = SPEAD_1D_daycounter(jday,iTime) 

jdaynew = ceil(iTime); % Current day
if jdaynew == jday
    newday = 'not';
else
    newday = 'yes';
    jday = jdaynew;
end

return
    
