function [B,stats,ychap] = myRobustFit(x,y)

[B,stats]=robustfit(x,y);

%........................
ychap = B(1) + B(2)*x;
%........................
% $$$ xmin=min(x(:));
% $$$ xmax=max(x(:));
% $$$ xchap = [xmin,(xmax-xmin)/100,xmax];
% $$$ ychap = B(1) + B(2)*xchap;
%........................
% $$$ figure(1)
% $$$ scatter(x,y,'filled','r')
% $$$ hold on
% $$$ plot(xchap,ychap,'b-')
% $$$ grid on
%........................

return
%**************************
%<http://www.mathworks.com/access/helpdesk/help/toolbox/stats/robustfit.html>

x = (1:10)';
y = 10 - 2*x + randn(10,1);
y(10) = 0;

bls = regress(y,[ones(10,1) x])

brob = robustfit(x,y)

scatter(x,y,'filled'); grid on; hold on
plot(x,bls(1)+bls(2)*x,'r','LineWidth',2);
plot(x,brob(1)+brob(2)*x,'g','LineWidth',2)
legend('Data','Ordinary Least Squares','Robust Regression')
