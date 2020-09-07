function [Apixel,AREAS]=mypixelarea() 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use: [Apixel,AREAS] = mypixelarea()
%
%Este programa calcula el Area [m2] de cada Pixel (1grad x 1grad) segun su latitud:
%OBTENGO LOS RADIOS (ri) DE LAS DIFERENTES CIRCUNFERENCIAS QUE APARECEN
%SI CORTO HORIZONTALMENTE (CON UN PLANO PARALELO AL ECUADOR) EL GLOBO
%TERRAQUEO A DISTINTAS LATITUDES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------
%Busco obtener el Arco (arcoi) que corresponde a 1grado de longitud segun me
%desplace en latitud. Ese arco va disminuyendo porque los meridianos se
%van juntando. Y en el fondo los meridianos se van juntando porque las
%circunferencias que se obtienen de cortar el Globo con un plano horizontal
%al ecuador para cada latitud son mas pequenas (menor ri).
%
%Entonces arcoi=f(ri). De modo que primero he obtener los ri que
%corresponden a cada latitud. Y segundo aplico la ecuacion del Arco para
%obtener cuantos metros mide 1grado de longitud a cada latitud.
%-------------------------------------------------------------------------
%DEF1: Obtencion del radio de las diferentes circunferencias que cortan
%el Globo a cada latitud:
%
%cos(alfa) = cateto contiguo/hipot ==> cateto contiguo = cos(alfa)*hipot
%
%Si tengo que:
%1. alfa = angulo de latitud (en radianes).
%2. hipot = radio de la Tierra (R)
%3. cateto contiguo = radio (ri) de la cicunferencia que corta el Globo Terrestre
%para cada latitud (ese radio va logicamente disminuyendo con la latitud
%hasta hacerse cero a 90grados; y es maximo a 0grados e igual al radio de
%la Tierra).
%
%Es decir:
%
%ri = cos(alfa)*R
%-------------------------------------------------------------------------
%DEF2: Ecuacion del Arco:
%
%360grados - 2*pi*ri (perimetro de la circunferencia)n
%1grado    -  arcoi (perimetro de un grado en la circunferencia)
%
%Es decir:
%
%1grado = (2*pi*ri)/360grados: conversion de 1grado a longitud (arco)
%segun el valor de ri (que varia con la latitud).
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PASO DE GRADOS A RADIANES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc=(2*pi)/360; %factor de conversion para pasar de Grados a Radianes: 1grado x (2*pi)rad/360grados = 0.017453 rad.
LAT=[0:90];
LATrad=fc*LAT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINO EL RADIO DE LA TIERRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=6378000 %[m] = 6.378 km (USAR EN m!!)
% $$$ R=6378 %[km]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO LOS ARCOS QUE CORRESPONDEN A 1grado DE LONG A CADA LATITUD:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RADIOS=[];
ARCOS=[];
for lati=LATrad
    %OBTENGO LOS RADIOS (ri) DE CADA CIRCUNFERENCIA HORIZONTAL:
    ri = cos(lati)*R;
    
    %OBTENGO LOS ARCOS QUE CORRESPONDEN A 1grado DE LONG EN CADA UNA DE
    %ESAS CIRCUNFERENCIAS HORIZONTALES (ES DECIR, A CADA LATITUD):
    arcoi = fc*ri;
    
    %STOCKO:
    RADIOS=[RADIOS;[lati/fc,ri]];
    ARCOS=[ARCOS,arcoi];
end
%RADIOS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMO LO QUE EN REALIDADA HE OBTENIDO SON LOS ARCOS A CADA LATITUD PERO
%YO TRABAJO CON PIXELS (ESTAN ENTRE DOS LATITUDES) HAYO LA MEDIA DE LOS
%ARCOS SUPERIOR E INFERIOR DE CADA PIXEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ainf=ARCOS(1:end-1);
Asup=ARCOS(2:end);
ARCpixel=(Ainf+Asup)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMO TENGO 2 HEMISFERIOS HAGO ESPEJO DE LO OBTENIDO:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ARC=[fliplr(ARCpixel),ARCpixel];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO EL AREA DE CADA PIXEL MULTIPLICANDO EL ARCO (lo puedo llamar la
%BASE) POR LA DISTANCIA DE 1grado EN LATITUD (lo puedo llamar la ALTURA)
%QUE ES INVARIANTE (NO CAMBIA SEGUN LA LATITUD):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
altura=fc*R; %arco (distancia) de 1grado de latitud.
AREAS=altura*ARC;
AREAS=AREAS(:); %(180,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONSTRUYO UNA MATRIZ (180,360):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Apixel=AREAS*ones(1,360); %(180,1) x (1,360) = (180,360): Area (en m^2) cubierta por cada pixel del Globo.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPRUEBO QUE ESTEN BIEN LAS AREAS DE CADA PIXEL (EL SUMATORIO DEBE
%DAR LA SUPERFICIE TOTAL DE LA TIERRA):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------
%La superficie de una esfera de radio, r, es   S = 4*pi*r^2
%------------------------------------------------------------
Aearth=4*pi*R^2;
sumApixel=sum(Apixel(:));
sumApixel/Aearth %debe salir igual a 1 para que sean correctas las areas de los pixeles calculadas.
return

%%%%%%%%
%PLOTEO:
%%%%%%%%
figure(1)
AApixel=Apixel;
AApixel(1,1)=-90;
imagesc(AApixel)
colorbar('horiz')

pause(1)
close(1)
