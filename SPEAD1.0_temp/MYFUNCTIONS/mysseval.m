function [sse] = mysseval(m,eqml,xdata,ydata)
%****************************************************************************
%Use [sse] = sseval_funany(m,eqml,xdata,ydata)
% m     -- vector of initial parameter values; that will go into eval(eqml)
% eqml  -- character string with the equation; eg. m(1)*exp(-((xdata-m(2))/m(3))^2)+m(4)
% xdata -- vector of independent variable observations.
% ydata -- vector of dependent variable observations f(x)
% ychap -- vector of dependent variable fitted curve fchap(x)
%****************************************************************************
ychap = eval(eqml);
sse_scale = 'linear';
%%sse_scale = 'logexp';
if     strcmp(sse_scale,'linear')
    ymean = mean(ydata);

    sst = sum(abs(ydata - ymean).^2); %total sum of squares (total variance VT)
    ssr = sum(abs(ychap - ymean).^2); %data sum of squares (explained variance VE)
    sse = sum(abs(ychap - ydata).^2); %error sum of squares (non-explained variance VNE)
    %%sst_bis = (sse + ssr); %(total variance VT)

% $$$     sse_bis = sum((ychap - ydata).^2); %Gives problems for negative parameters (ie. ychap has imaginary numbers) -- use the "abs"
% $$$     if abs(sse - sse_bis) > 1d-3
% $$$ 	m
% $$$ 	xdata
% $$$ 	ychap
% $$$ 	ydata
% $$$ 	ydist = (ychap - ydata).^2
% $$$ 	sse
% $$$ 	sse_bis
% $$$ 	pause
% $$$     end
    
elseif strcmp(sse_scale,'logexp')

    ychap_log = log(ychap);
    ydata_log = log(ydata);
    ymean_log = mean(ydata_log);

    sst = sum(abs(ydata_log - ymean_log).^2); %total sum of squares (total variance VT)
    ssr = sum(abs(ychap_log - ymean_log).^2); %data sum of squares (explained variance VE)
    sse = sum(abs(ychap_log - ydata_log).^2); %error sum of squares (non-explained variance VNE)
    %%sst_bis = (sse + ssr); %(total variance VT)

end

rsquare = 1.0 - (sse/sst);
%%rsquare_bis = ssr/(sse+ssr);
rhocorr = sqrt(rsquare);

% $$$ sst
% $$$ ssr
% $$$ sse
% $$$ sst_bis %SHOULD BE THE SAME BUT THEY ARE NOT!!!
% $$$ rsquare
% $$$ rsquare_bis %SHOULD BE THE SAME BUT THEY ARE NOT!!!
% $$$ rhocorr
% $$$ pause

return
