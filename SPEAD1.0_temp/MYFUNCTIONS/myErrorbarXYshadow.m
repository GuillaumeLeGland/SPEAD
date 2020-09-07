function [hpatch]=myErrorbarXYshadow(x,y,stdx,stdy)

%=======================================
%..............
stdx=stdx(:);
stdy=stdy(:);
nptos = length(stdx);
%..............
hpatch=[];    
for i=1:nptos
    xi=x(i);
    yi=y(i);
    stdxi=stdx(i);
    stdyi=stdy(i);
    %......
    xy=[xi,yi]
    stdxy=[stdxi,stdyi]
    %......
    hpatchi=myElipticPatch(xi,yi,stdxi,stdyi);
    hpatch=[hpatch,hpatchi];
    pause
end
% $$$ hold on
%..............
% $$$ hplot = plot(x,y,'*');
% $$$ hold off
%..............
% $$$ set(gca,'Xlim',[0 1])
% $$$ set(gca,'Ylim',[0 1])
%..............
return
%=======================================
function [CH]=myElipticPatch(x0,y0,Xradius,Yradius)

t=0:0.01:2*pi;
x=x0+Xradius*cos(t);
y=y0+Yradius*sin(t);

%.....................
% $$$ myColor=[0.5 0.5 0.5];
% $$$ myTransparency=0.2;
%.....................
myColor=[0.9 0.9 0.9];
myTransparency=1.0;
%.....................

CH = patch(x,y,myColor);
set(CH,'facealpha',myTransparency,'EdgeColor','none');

return
%=======================================

%********************************
return

%EXAMPLE:
X=rand(3,4,100);
Y=rand(3,4,100);

meanX=mean(X,3);
meanY=mean(Y,3);

stdX=std(X,[],3)/10;
stdY=std(Y,[],3)/10;

x=meanX(:);
y=meanY(:);

stdx=stdX(:);
stdy=stdY(:);

figure(1)
subplot(2,2,1)
myErrorbarXYshadow(meanX,meanY,stdX,stdY)
subplot(2,2,2)
myErrorbarXYshadow(x,y,stdx,stdy)

