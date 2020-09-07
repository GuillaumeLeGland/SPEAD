function []=mymonthlymeansbyoceanregionplots(VARNAME,fignum,varargin) %eq.to:codim_byregionplots(VARNAME,TERRA,MAPterra,fignum,varargin).

%**************************************
%Use: []=mymonthlymeansbyoceanregionplots(VARNAME,fignum,VAR1,VAR2,...etc)
%**************************************
n=length(varargin);

[ZONAS,TERRA,MAPterra]=myoceanregions;

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

    mVARSzonai=[];
    for vari=1:n
	mVAR=varargin{vari};
	mVARzonai=mVAR(i,:);
	mVARSzonai=[mVARSzonai;mVARzonai];
    end
    
    %æææææææææææææææææææææææææææææææææææææææææææææææ
    %PLOTS CCN vs CCNbio MEDIAS TOTALES (12 ptos):
    %æææææææææææææææææææææææææææææææææææææææææææææææ
    figure(fignum)
    subplot(4,4,i)
    if strcmp(VARNAME,'DMS')
	ymin=0;
	ymax=12;
	dy=2;
	mVARzonai
	hp=plot([1:12],mVARSzonai(1,:),'-',[1:12],mVARSzonai(2,:),'-',...
		[1:12],mVARSzonai(3,:),'-',[1:12],mVARSzonai(4,:),'-',...
		[1:12],mVARSzonai(5,:),'-',[1:12],mVARSzonai(6,:),'-',...
		[1:12],mVARSzonai(7,:),'-');
% $$$ 		[1:12],mVARSzonai(7,:),'-',[1:12],mVARSzonai(8,:),'-');
	set(hp(1),'Color','red','Marker','+','Markersize',[5]);
	set(hp(2),'Color',[0 0.8 0],'Marker','x','Markersize',[5]);
	set(hp(3),'Color','blue','Marker','.','Markersize',[5]);
	set(hp(4),'Color','cyan','Marker','*','Markersize',[5]);
	set(hp(5),'Color','magenta','Marker','o','Markersize',[3]);
	set(hp(6),'Color',[0.8 0 0],'Marker','d','Markersize',[5]);
	set(hp(7),'Color',[0 0 0],'Marker','>','Markersize',[5]);
% $$$ 	set(hp(8),'Color',[0.5 0.5 0.5],'Marker','<','Markersize',[5]);
    elseif strcmp(VARNAME,'DMSF')
	ymin=0;
	ymax=35;
	dy=5;
	mVARSzonai
	hp=plot([1:12],mVARSzonai(1,:),'-',[1:12],mVARSzonai(2,:),'-',...
		[1:12],mVARSzonai(3,:),'-',[1:12],mVARSzonai(4,:),'-',...
		[1:12],mVARSzonai(5,:),'-');
% $$$ 		[1:12],mVARSzonai(5,:),'-',[1:12],mVARSzonai(6,:),'-');
	set(hp(1),'Color','blue','Marker','.','Markersize',[5]);
	set(hp(2),'Color','cyan','Marker','*','Markersize',[5]);
	set(hp(3),'Color','magenta','Marker','o','Markersize',[3]);
	set(hp(4),'Color',[0.8 0 0],'Marker','d','Markersize',[5]);
	set(hp(5),'Color',[0 0 0],'Marker','>','Markersize',[5]);
% $$$ 	set(hp(6),'Color',[0.5 0.5 0.5],'Marker','<','Markersize',[5]);
    elseif strcmp(VARNAME,'CHL')
	ymin=0;
	ymax=2.5;
	dy=0.5;
	hp=plot([1:12],mVARSzonai(1,:),'-',[1:12],mVARSzonai(2,:),'-',...
		[1:12],mVARSzonai(3,:),'-',[1:12],mVARSzonai(4,:),'-',...
 		[1:12],mVARSzonai(5,:),'-',[1:12],mVARSzonai(6,:),'-');
	set(hp(1),'Color','red','Marker','+','Markersize',[5]);
	set(hp(2),'Color','blue','Marker','.','Markersize',[5]);
	set(hp(3),'Color','cyan','Marker','*','Markersize',[5]);
	set(hp(4),'Color','magenta','Marker','o','Markersize',[3]);
	set(hp(5),'Color',[0.8 0 0],'Marker','d','Markersize',[5]);
 	set(hp(6),'Color',[0 0 0],'Marker','>','Markersize',[5]);
    end
    posy=ymin+(0.03*(ymax-ymin));
    posx=2;
    %.................................................................................
    set(gca,'Xlim',[0 13],'XTick',[1:1:12],'XTickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'])
    set(gca,'Ylim',[0 ymax],'YTick',[0:dy:ymax],'Ycolor',[0 0 0]);
    set(gca,'FontSize',[4]);
    grid on
    %.................................................................................
    ht=title([num2str(i),': ',zona]);
    set(ht,'Fontsize',[6])
    if i==15
	if strcmp(VARNAME,'DMS')
	    hl=legend('Kettle','Simo','POP','PISCES','PLANKTOM','HAMOCC','HADOCC');
% $$$ 	    hl=legend('Kettle','Simo','POP','PISCES','PLANKTOM','HAMOCC','HADOCCander','HADOCCsimo');
	    %get(hl,'Position')
	    set(hl,'Position',[0.6,0.19666,0.057031,0.06857])
	elseif strcmp(VARNAME,'DMSF')
	    hl=legend('POP','PISCES','PLANKTOM','HAMOCC','HADOCC');
% $$$ 	    hl=legend('POP','PISCES','PLANKTOM','HAMOCC','HADOCCander','HADOCCsimo');
	    %get(hl,'Position')
	    set(hl,'Position',[0.6,0.19666,0.057031,0.06857])
	elseif strcmp(VARNAME,'CHL')
	    hl=legend('SeaWiFS','POP','PISCES','PLANKTOM','HAMOCC','HADOCC');
	    %get(hl,'Position')
	    set(hl,'Position',[0.57,0.19666,0.057031,0.06857])
	end
	set(hl,'Fontsize',[6])
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
