function [ychap,r2] = mylinearegression(xdata,ydata)
    
x = xdata(:); 
y = ydata(:); 

xr = x; 

X  = [ones(length(x),1),x];
Xr = [ones(length(xr),1),xr];

B = (X'*X)\(X'*y); %coefficientes de regression.

ychap = Xr*B; %recta linear de regression.

ym = sum(y)/length(y);

vt  = sum((y-ym).^2); 
ve  = sum((ychap-ym).^2); 
vne = sum((y-ychap).^2); 

r2 = ve/vt; %Coefficient of determination. 

return

