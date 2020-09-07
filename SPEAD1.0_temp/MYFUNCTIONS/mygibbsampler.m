%***********************************************
% <http://timsalimans.com/the-power-of-jit-compilation/>
% Use:
% nptos = 5e4;
% thin = 1e3;
% [xs,ys] = mygibbsampler(nptos,thin);
% figure(1)
% plot(xs,ys,'.')
%***********************************************

function [x_samp,y_samp] = mygibbsampler(nptos,thin)
feature('accel','on') %JIT accelarator -- disabled (off) or enabled (on)
x_samp = zeros(nptos,1);
y_samp = zeros(nptos,1);
y=0;
for i=1:nptos
    gammarands=randg(3,thin,1);
    normrands=randn(thin,1);
    for j=1:thin
        x=(y^2+4)*gammarands(j);
        y=1/(1+x)+normrands(j)/sqrt(2*x+2);
    end
    x_samp(i) = x;
    y_samp(i) = y;
end
return

