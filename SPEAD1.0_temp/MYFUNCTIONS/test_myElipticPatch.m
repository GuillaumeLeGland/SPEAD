x0=5;
y0=5;

Xradius=0.50;
Yradius=0.25;

t=0:0.01:2*pi;
x=x0+Xradius*cos(t);
y=y0+Yradius*sin(t);

myColor='b';
myTransparency=0.2;

CH = patch(x,y,myColor);
set(CH,'facealpha',myTransparency,'edgecolor',myColor,'EdgeColor','none');
set(gca,'Xlim',[0 10])
set(gca,'Ylim',[0 10])
