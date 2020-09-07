function [mX]=mymonthlymeansbyoceanregion(X) %equivalent to '[mX]=codim_byregion(X,ZONAS)'.

%**************************************************
%Use: [mX]=mymonthlymeansbyoceanregion(X,ZONAS).
%
%Need first to use: [ZONAS,TERRA,MAPterra]=myoceanregions().
%**************************************************

[ZONAS,TERRA,MAPterra]=myoceanregions;
    
%%%%%%%%%%%%%%%
%MATRIZ SIZE:
%%%%%%%%%%%%%%%
[m,n,p]=size(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QUITO POSIBLES VALORES NEGATIVOS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=find(X<0);
X(I)=nan;

%%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO MONTHLY MEANS:
%%%%%%%%%%%%%%%%%%%%%%%
mX=[];
for i=1:15
    zona=i;
    zonai=ZONAS(i,:); %posiciones en el WordMap de cada zonai.
    I=find(isfinite(zonai));
    zonai=zonai(I); %quito los NaN.
    if i==9
	zonaMedit=zonai; %(la usare para el mapa con las estaciones)
    end

    %===============================
    %MONTHLY-MEANS PARA CADA REGION:
    %===============================
    mXzonai=[];
    for k=1:p
	Xk=X(:,:,k);
	%VALORES (PIXELS) PARA LA REGION:
	Xzonaik=Xk(zonai);
	%VALORES MEDIOS PARA LA REGION:
	mXzonaik=nanmean(Xzonaik);
	%STOCKAGE:
	mXzonai=[mXzonai;mXzonaik]; %medias mensuales para la zonai
    end%endif k=1:p

    %=========
    %STOCKAGE:
    %=========
    mX=[mX;mXzonai(:)']; %(15zonas x 12meses)
end %end esti=1:15

