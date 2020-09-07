function [r,t,varargout]=mycorrcoef(X,Y,Mode)
%*********************************************************************
%PROGRAMA "CORR_CORRCOEF.m": Este programa calcula el coeficiente de
%correlacion (PEARSON o SPEARMAN, segun se indique en "Mode") entre 2
%vectores X e Y. No importa que haya NaN's.
%    
%Use: [r,t]=mycorrcoef(X,Y,'mode').
%   
%donde:
% X e Y: son los vectores a correlacionar.
% 'mode': es un string que indique el metodo de corr: 'pearson' o 'spearman'.
% r: coef. de correlacion.
% t: significancia (debe ser inferior a -1.9 o superior a +1.9 para ser un valor de corr. significativo)
% p: p-value.
%*********************************************************************
%----------------------------------------------------------------------
%NOTA: Hacer la correlacion de SPEARMAN, en realidad no es mas que
%aplicar la formula de la correlacion de PEARSON al ranking (posicion
%relativa en el vector segun su valor) de mis datos, en lugar de a los
%datos (el valor) en si. De esta manera, se disminuye mucho el problema
%de los outlayers (ya no se juega con valores sino con posiciones
%relativas de esos valores unos respecto a otros, en el vector) y NO se
%hacen ningun tipo de asumpcion parametrica (ej. Gaussianidad,
%bigaussianidad, homogeneidad de varianzas, etc).
%----------------------------------------------------------------------    

%PONGO LOS VECTORES EN COLUMNA:
X=X(:);
Y=Y(:);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PONGO EN NaN EN LAS POSICIONES DE LOS DOS VECTORES DONDE HAYA NaN PARA
%UNO DE ELLOS (SE NECESITAN PAREJAS DE DATOS PARA CALCULAR LA CORRELACION):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inan=find(isnan(X) | isnan(Y));
X(Inan)=nan;
Y(Inan)=nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPRUEBO QUE HAY PAREJAS DE DATOS Y SI ES ASI CALCULO LA CORRELACION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=find(isnan(X)==0 & isinf(X)==0); %donde NO hay NaN's o Inf's.
J=find(isnan(Y)==0 & isinf(X)==0);

if length(I)>0 & I==J %tengo datos dentro de los 2 vectores (en las misma pos).
    X=X(I);
    Y=Y(I);

    if strcmp(lower(Mode(1:7)),'pearson') %PEARSON

        if length(X)==length(Y), N=length(X);
        else display('X e Y deben tener igual longitud'),pause
        end
	
        %ææææææææææææææææææææææææææææææææææ
        %CALCULO EL COEF DE CORR. PEARSON:
        %ææææææææææææææææææææææææææææææææææ
        covxy=(1/(N-1))*(sum((X-mean(X)).*(Y-mean(Y))));
        r=covxy/(std(X)*std(Y));
	
    elseif strcmp(lower(Mode(1:8)),'spearman'); %SPEARMAN
	%ææææææææææææææææææææææææ
        %ORDENO DE MENOR A MAYOR:
	%ææææææææææææææææææææææææ
        [v1,e1]=sort(X); %e1=etiqueta (pos.original de los valores en el vect X)
        [v2,e2]=sort(Y);
        [x1,RG1]=sort(e1); %obtengo los rangos (ranks) (la x1 no se usa)
        [x2,RG2]=sort(e2);
        N=length(v1); %meses con dato.

	%æææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
        %COMPRUEBO QUE NO HAYA NINGUN rg REPETIDO (DEBIDO A QUE SE REPITA UN
        %VALOR EN EL VECTOR X o Y):
	%æææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
        %-------------------------------------------------------------
        %Si un mismo valor se encuentra en varias posic (esta repetido), 
        %asigno el rg medio (ver pag. 252 del W.J. CONOVER, "Practical
        %Non-Parametric Statistics")
        %-------------------------------------------------------------
	diffv1=diff(v1);
	diffv2=diff(v2);
	
        I1=find(diffv1==0);
        if length(I1)>=1
	    %HAGO LAS MEDIAS POR GRUPOS DE DATOS REPETIDOS (PUEDE HABER
            %MAS DE UN VALOR QUE SE REPITA):
	    diffv1=[diffv1;nan]; %anado el nan para el "if diffv1(i)==diffv1(i+1)"
	    for i=I1'
		%SI HAY 3 O MAS DATOS REPES SEGUIDOS (diffv1 TIENE 2 O MAS CEROS CONTIGUOS) 
		%COJO EL ULTIMO (ES PARA NO MALGASTAR CALCULOS):
		if diffv1(i)==diffv1(i+1)
		    continue
		end
		H1=find(v1==v1(i)); %pos en v1 donde esta un grupo de datos repes.
		J1=e1(H1); %obtengo las etiquetas que corresponden.
		rgm1=sum(RG1(J1))/length(J1); %hayo el rango medio.
		RG1(J1)=rgm1;
	    end
        end
        I2=find(diffv2==0);
        if length(I2)>=1
	    diffv2=[diffv2;nan];
	    for i=I2'
		if diffv2(i)==diffv2(i+1)
		    continue
		end
	        H2=find(v2==v2(i)); 
		J2=e2(H2);
		rgm2=sum(RG2(J2))/length(J2);
		RG2(J2)=rgm2;
	    end
        end

% $$$ 	whos RG1 RG2
% $$$ 	[X(:),RG1(:),Y(:),RG2(:)]
% $$$ 	pause
        
	%ææææææææææææææææææææææææææææææææææ
        %CALCULO EL COEF DE CORR. SPEARMAN:
        %ææææææææææææææææææææææææææææææææææ
        %d=RG1-RG2; %dif. entre las posic que ocupan los valores de cada estimador
        %d2=d.^2;     
        %r=1-6*(sum(d2)/(N*((N^2)-1))); %correlacion spearman
        r = 1 - ( 6*sum((RG1-RG2).^2) / (N*(N^2-1)) ); %correlacion spearman
 	rbis = 1-6*sum((RG1-RG2).^2)/N/(N^2-1);

	%ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
        %CALCULO EL COEF DE CORR. PEARSON (PERO SOBRE LOS RANKINS):
        %ææææææææææææææææææææææææææææææææææææææææææææææææææææææææææ
        covxy=(1/(N-1))*(sum((RG1-mean(RG1)).*(RG2-mean(RG2))));
        rr=covxy/(std(RG1)*std(RG2));
	
% $$$   	[N,r,rr,rbis],pause(1)
	%.......................................................................
	%test: OJO! ESTE TEST NO SE CUMPLE SIEMPRE CUANDO QUIERO CALCULAR LA CORR "GLOBAL"
        %DE SPEARMAN SOBRE DMS Y CCN (corr16.m). NO SE PQ, COMPROBARLO...
	%(CREO QUE ES ALGO RELACIONADO CON EL NUMERO DE ZEROS EN LA
        %MATRIZ... MEJOR USAR MATRICES CON NaN)
	deltarho=r-rr;
	if deltarho>10^(-3)
	    deltarho
	    [r,rr]
	    display('Error: r y rr provienen de formulas equivalentes, deben ser iguales')    
	end
	%........................................................................
    end %endif strcmp
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SIGNIFICANCIA: (Test "t": Ho: nu=0 ; H1: nu~=0;)
    %Sr = sqrt((1-r^2)/(n-2)); 
    %t = (r - nu)/Sr = (r - 0)/Sr = r*sqrt((n-2)/(1-r^2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if abs(r)~=1
        t = r*sqrt((N-2)/(1-r^2));
    elseif abs(r)==1
        t = inf; %asi me quito el 'Warning Divided by Zero'
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the p-values:
    %%%%%%%%%%%%%%%%%%%%%%%%%
    p=2*(1-mytcdf(abs(t),N-2));
    
else %todo el vector son NaN's
    r=nan;
    t=nan;
    p=nan;
end
nrtp=[N,r,t,p];
varargout{1}=N;
varargout{2}=p;
