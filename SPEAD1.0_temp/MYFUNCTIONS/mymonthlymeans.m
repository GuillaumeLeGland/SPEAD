function [mX,stdX,Nptos]=mymonthlymeans(x,time)

%***************************************************************************
%Program MYMONNTHLYMEANS.m: Este programa calcula las monthly means a
%partir de un vector "x" con datos para todo un ano (del dia 1 al dia 365).
%
%Use: [mX,stdX,nptos]=mymonthlymeans(x,time)
%
%donde:
%
% x: Vector original de datos para todo un ano.
% time: Vector con los dias al que corresponde cada dato.
% mX: 12 monthly means.
% stdX: 12 standart deviations de las monthly means.
%***************************************************************************

%..................................    
% $$$ monthlims=[0,31,59,90,120,151,181,212,243,273,304,334,365];
% $$$ monthdays=[31,28,31,30,31,30,31,31,30,31,30,31];
%..................................    
monthlims=[0:30:360];
monthdays=[30,30,30,30,30,30,30,30,30,30,30,30];
%..................................    
S=size(x); %[lat,lon,time]
H=find(S>1);
ndim=length(H);
%..................................    
if ndim==1 %1D vector [time]
    %..............................
    m=length(x);
    n=length(time);
    if m~=n
	display('error, "x" y "time" deben ser de igual length!')
    end
    %..............................
    mX=[];
    stdX=[];
    Nptos=[];
    for k=1:12
	mesk=k;
	day1=monthlims(k);
	day30=monthlims(k+1);
	I=find(time>day1 & time<=day30); %aprox. 30 days con dato para el mesk.
	xk=x(I); %datos para el mesk.
	H=find(isnan(xk)==0); %donde hay datos.
	nptos=length(H);
	if nptos>=1 %meses con datos.
	    mx=nanmean(xk);
	    stdx=nanstd(xk);
	else %meses sin ningun dato.
	    mx=nan;
	    stdx=nan;
	    nptos=nan;
	end
	mX=[mX,mx];
	stdX=[stdX,stdx];
	Nptos=[Nptos,nptos];
    end
    %........................
% $$$     figure(200)
% $$$     months=[15:30:365]; %dia en el medio de cada mes.
% $$$     plot(time,x,'-b.')
% $$$     hold on
% $$$     plot(months,mX,'-r.')
% $$$     pause(1)
% $$$     close(200)
    %........................
elseif ndim==3 %3D array [lat,lon,time]
    %..............................
    [mm,nn,pp]=size(x);
    p=length(time);
    if pp~=p
	display('error, "x" y "time" deben ser de igual length!')
    end
    %..............................
    mX=[];
    stdX=[];
    Nptos=[];
    for k=1:12
	mesk=k;
	%..............
	day1=monthlims(k);
	day30=monthlims(k+1);
	%..............
	I=find(time>day1 & time<=day30); %aprox. 30 days con dato para el mesk.
	xk=x(:,:,I); %datos para el mesk.
	%..............
	mx=nanmean(xk,3);
	stdx=nanstd(xk,0,3);
	%..............
	mX(:,:,k)=mx;
	sdtX(:,:,k)=stdx;
	%..............
	Nptos=nan;
	%..............
    end
    %..............................
% $$$     figure(200)
% $$$     for k=1:12
% $$$ 	Xk=mX(:,:,k);
% $$$ 	subplot(3,4,k)
% $$$ 	imagesc(Xk)
% $$$ 	colorbar
% $$$     end
    %..............................
end

