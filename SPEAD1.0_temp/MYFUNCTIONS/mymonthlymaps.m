function []=mymonthlymaps(VAR,VARNAME,Land,palete,fignum)

%************************************
%Programa MYMONTHLYMAPS.m:
%Este programa grafica 12 subplots con monthly-means values.
%Uso: []=mymonthlymaps(VAR,'VARNAME',Land,palete,fignum)
%************************************

VARNAME %tiene que ser un 'string'
meses=['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
%...............................................................
gradoslat=[+90,+80,+60,+40,+20,0,-20,-40,-60,-80,-90];
gradoslong=[-180,-160,-140,-120,-100,-80,-60,-40,-20,0,+20,+40,+60,+80,+100,+120,+140,+160,+180];
lattick=[1,[10:20:170],180];
longtick=[1,[20:20:360]];
%...............................................................
gradoslat2=[+80,+60,+40,+20,0,-20,-40,-60,-80];
gradoslong2=[-160,-120,-80,-40,0,+40,+80,+120,+160];
lattick2=[10:20:170];
longtick2=[20:40:340];
%...............................................................

maxVAR=max(VAR(:));
minVAR=min(VAR(:));
%pcnt=0.05;
pcnt=0.025;
limsup=maxVAR+(pcnt*maxVAR);
if minVAR>=0
    liminf=minVAR-(pcnt*maxVAR);
elseif minVAR<0
    liminf=minVAR-(pcnt*abs(minVAR));
    liminf=-max([abs(limsup),abs(liminf)]);
end
limsup,liminf

figure(fignum)
MAP=palete;
MAP=[[0 0 0];MAP;[0.5 0.5 0.5]];
for k=1:12
    subplot(3,4,k)
    mesk=meses(k,:);
    VARk=VAR(:,:,k);
    VARk(Land)=limsup;
    VARk(1,1)=liminf;
    VARk(1,2)=limsup;
    imagesc(VARk);
    hc=colorbar('horiz');
    colormap(MAP)
    ht=title([VARNAME,' ',mesk]);
% $$$     set(gca,'YTick',lattick,'YTickLabel',gradoslat,'XTick',longtick,'XTickLabel',gradoslong,'FontSize',[6])
    set(gca,'YTick',lattick2,'YTickLabel',gradoslat2,'XTick',longtick2,'XTickLabel',gradoslong2,'FontSize',[4])
    set(hc,'FontSize',[6]);
    set(ht,'FontSize',[6])
    grid on
end

