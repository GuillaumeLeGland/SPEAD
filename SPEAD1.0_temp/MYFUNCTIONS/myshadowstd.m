x=[1:12];
y=rand(12,1);
stdy=0.1*y;
z1=y+stdy;
z2=y-stdy;

figure(1)
a1=area(x,z1);
grid on
set(gca,'Layer','top')
set(a1,'LineStyle','none');
set(a1,'FaceColor',[0.9 0.9 0.9]);
hold on;
a2=area(x,z2);
grid on
set(gca,'Layer','top')
set(a2,'LineStyle','none');
set(a2,'FaceColor',[1 1 1]);
hold on;
