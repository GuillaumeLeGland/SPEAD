function []=mymonthlymeansbyoceanregionplotsXYZ(VARNAMEx,VARNAMEy,VARNAMEz,fignum,X,Y,Z,xscalemax,dx) 

%**************************************
%Use: []=mymonthlymeansbyoceanregionplotsXYZ('VARNAME1','VARNAME2','VARNAME3',fignum,VAR1(180,360,12),VAR2(180,360,12),VAR2(180,360,12),xscalemax,dx)
%**************************************

[ZONAS,TERRA,MAPterra]=myoceanregions;

[mX]=mymonthlymeansbyoceanregion(X); %(15,12)
[mY]=mymonthlymeansbyoceanregion(Y);
[mZ]=mymonthlymeansbyoceanregion(Z);

%====================================================
%OBTENGO LA CORRELACION DE SPEARMAN PARA CADA REGION ENTRE X e Y (dejo fuera Z):
%====================================================
RS=[];
for i=1:15
    mXzonai=mX(i,:);
    mYzonai=mY(i,:);

    [rsx,tsx]=mycorrcoef(mXzonai,mYzonai,'Spearman');
    RS=[RS;rsx];
end%endif i=1:15 zonas.

%=========
%GRAFICAS:
%=========
for i=1:15
    if i==1
	zona='North Atlantic (40^{\circ}N - 70^{\circ}N)';
    elseif i==2
	zona='Subtrop. Atlantic NH (10^{\circ}N - 40^{\circ}N)';
    elseif i==3
	zona='Equat. Atlantic (10^{\circ}S - 10^{\circ}N)';
    elseif i==4
	zona='Subtrop. Atlantic SH (10^{\circ}S - 40^{\circ}S)';
    elseif i==5
	zona='North Pacific (40^{\circ}N - 70^{\circ}N)';
    elseif i==6
	zona='Subtrop. Pacific NH (10^{\circ}N - 40^{\circ}N)';
    elseif i==7
	zona='Equat. Pacific (10^{\circ}S - 10^{\circ}N)';
    elseif i==8
	zona='Subtrop. Pacific SH (10^{\circ}S - 40^{\circ}S)';
    elseif i==9
	zona='Mediterranean Sea';
    elseif i==10
	zona='Subtrop. Indian NH (10^{\circ}N - 40^{\circ}N)';
    elseif i==11
	zona='Equat. Indian (10^{\circ}S - 10^{\circ}N)';
    elseif i==12
	zona='Subtrop. Indian SH (10^{\circ}S - 40^{\circ}S)';
    elseif i==13
	zona='Arctic (70^{\circ}N - 90^{\circ}N)';
    elseif i==14
	zona='Southern Ocean (40^{\circ}S - 60^{\circ}S)';
    elseif i==15
	zona='Antarctic (60^{\circ}S - 90^{\circ}S)';
    end
    mXzonai=mX(i,:);
    mYzonai=mY(i,:);
    mZzonai=mZ(i,:);

    %æææææææææææææææææææææææææææææææææææææææææææææææ
    %PLOTS CCN vs CCNbio MEDIAS TOTALES (12 ptos):
    %æææææææææææææææææææææææææææææææææææææææææææææææ
    xmin=0;
    ymin=0;
    xmax=max(mX(:));
    ymax=max(mY(:));
    posx=xmin+(0.03*(xmax-xmin));
    posy=ymin+(0.03*(ymax-ymin));

    figure(fignum)
    meses=['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
    subplot(4,4,i)
    hp=plot([1:12],mXzonai,'xb--',[1:12],mYzonai,'.r-',[1:12],mZzonai,'*k');
    set(gca,'Xlim',[0 13],'XTick',[1:1:12],'XTickLabel',meses,'FontSize',[5]);
    set(gca,'Ylim',[0 xscalemax],'YTick',[0:dx:xscalemax],'YTickLabel',[0:dx:xscalemax]);
% $$$     set(hp(3),'Color',[0 0.5 0])
    ht=title([num2str(i),': ',zona]);
    set(ht,'Fontsize',[6])
    grid on
    if i==4 %i==15
	hl=legend([hp(1),hp(2),hp(3)],{VARNAMEx,VARNAMEy,VARNAMEz});
	set(hl,'Fontsize',[5])
    end
    
    %ANADO VALOR DE CORRELACION:
    rsx=RS(i);
    corr=num2str(round(100*rsx)/100);
    %htxt=text(1.1,0.8*xscalemax,['\rho = ',corr]);
    htxt=text(0.25,0.8*xscalemax,['\rho = ',corr]);
    set(htxt,'Fontsize',[8],'Color',[0 0 0])

    if i==15
	%ANADO MAPA DE COLORES DE REGIONES:
	subplot(4,4,16)
	imagesc(TERRA)
	hc=colorbar('horiz');
	colormap(MAPterra)
	imagesc(TERRA) %repito para que se coloque bien el colormap.
	hc=colorbar('horiz');
	colormap(MAPterra)
	set(hc,'Position',[0.748845 0.1 0.156155 0.0123161],'XTick',[0.5:14.5],'XTickLabel',[1:15],'Fontsize',[4])
	gradoslat2=[+80,+60,+40,+20,0,-20,-40,-60,-80];
	gradoslong2=[-160,-120,-80,-40,0,+40,+80,+120,+160];
	set(gca,'YTick',[10:20:170],'YTickLabel',gradoslat2,'XTick',[20:40:340],'XTickLabel',gradoslong2,'FontSize',[4])
    end
end

