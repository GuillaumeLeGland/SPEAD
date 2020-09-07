function [a,b,bs,r2,ychap,nptos] = myRegression(y_vector,X_vectors) %[corte-eje,pente,coef-determinacion]
%****************************************************************************
%Use: [a,b,bs,r2,ychap,nptos] = myRegression(y_vector,X_vectors) %[corte-eje,pente,coef-determinacion]
%****************************************************************************

y = y_vector(:);
numvect=min(size(X_vectors)); %numero de vectores.    

%PASO A VECTORES COLUMNA:
[j]=find(size(X_vectors)==min(size(X_vectors)));
if j==1
    X_vectors=X_vectors';
end

%COMPRUEBO QUE SON VECTORES DE IGUAL LONGITUD:
my=length(y);
mx=size(X_vectors,1);
if mx~=my, error('Los vectores "x1,x2..." e "y" deben ser de igual longitud')
end

if numvect==1 %REGRESION SIMPLE (2 VARIABLES)    
    x=X_vectors;
 
    %COMPRUEBO QUE SON VECTORES DE IGUAL LONGITUD:
    mx=length(x);
    if mx~=my, error('Los vectores x e y deben ser de igual longitud')
    end

    %QUITO NaN:
    I=find(isfinite(x) & isfinite(y));
    x=x(I);
    y=y(I);
    nptos=length(x);

    %COEFS. DE LA REGRESION:
    x=x(:);%pongo en columna.
    y=y(:);
    X=[ones(length(x),1),x];
    B=(X'*X)\(X'*y);

    %CORTE EN EL EJE Y PENDIENTE:
    a=B(1);
    b=B(2);

    %PENDIENTE DE LA REGRESION STANDARIZADA (COINCIDE CON EL COEFICIENTE DE
    %CORRELACION EN LAS REGRES SIMPLE).
    bs=b*std(x)/std(y);
    
    %RECTA DE REGRESION:
    ychap=X*B;

    %R2 (COEF. DETERMINACION):
% $$$     ym=sum(y)/length(y);
    ym=mean(y);
    vt=sum((y-ym).^2);
    ve=sum((ychap-ym).^2);
    vne=sum((y-ychap).^2);
    r2=ve/vt;

    %PONGO LOS NaN DONDE ESTABAN:
    yychap=ones(my,1)*nan;
    yychap(I)=ychap;
    ychap=yychap;

    %PLOTEO:
    %%plot(x,y,'r.',x,ychap,'b-')

elseif numvect>1

    %QUITO NaN:
    [I,J]=find(isnan([y,X_vectors]));
    X_vectors(I,:)=[];
    y(I)=[];
    nptos=length(X_vectors(:,1));

    %COEFS. DE LA REGRESION:
    X=[ones(nptos,1),X_vectors];
    B=(X'*X)\(X'*y);
    
    %CORTE EN EL EJE Y PENDIENTES:
    a=B(1);
    b=B(2:end);

    %COEFS. STANDARIZADOS DE LA REGRESION (NO! CORRESPONDEN CON LOS
    %COEFICIENTES DE CORRELACION EN LAS REGRES MULTIPLES!! OJO!)
    bs=b.*(std(X_vectors)/std(y))';
    
    %RECTA DE REGRESION:
    ychap=X*B;
    
    %R2 (COEF. DETERMINACION):
    ym=mean(y);
    ve = norm(ychap-ym)^2;  %Varianza explicada.
    vt = norm(y-ym)^2;      %Varianza total.
    vne= norm((y-ychap).^2);%Varianza no-explicada.
    r2=ve/vt;

    %PONGO LOS NaN DONDE ESTABAN:
    J=find(isfinite(y_vector));
    yychap=ones(my,1)*nan;
    yychap(J)=ychap;
    ychap=yychap;

    %................................................................
    %VISUALIZON EN 3D (SOLO SI HAY 2 VARIABLES INDEPENDIENTES):
% $$$     if numvect==2 %Y depende de 2 variables (Y = b0 + b1X1 + b2X2)
% $$$ 	[X1,X2]=meshgrid(X_vectors(:,1),X_vectors(:,2));
% $$$ 	Ychap = b(1) + b(2).*X1 + b(3).*X2;
% $$$ 	
% $$$ 	figure(1)
% $$$ 	surf(X1,X2,Ychap)
% $$$ 	pause
% $$$     end
    %................................................................
end
%****************************
return
%%disp('End of "myRegress.m"')
