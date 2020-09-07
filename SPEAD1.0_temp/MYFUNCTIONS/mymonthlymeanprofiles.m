function [mX,stdX]=mymonthlymeanprofiles(X)

%***************************************************************************
%Program MYMONNTHLYMEANS.m: Este programa calcula los monthly mean profiles a
%partir de una matriz 2D(m,365) "X" con datos para todo un ano (del dia 1 al dia 365).
%
%Use: [mX,stdX,nptos]=mymonthlymeanprofiles(X)
%
%donde:
%
% X: Matriz 2D(depth,time) con datos para todo un ano.
% time: Vector con los dias al que corresponde cada dato.
% mX: 12 monthly means.
% stdX: 12 standart deviations de las monthly means.
%***************************************************************************
[m,n]=size(X);
monthlims=[1,31,59,90,120,151,181,212,243,273,304,334,365];
monthdays=[31,28,31,30,31,30,31,31,30,31,30,31];
mX=[];
stdX=[];
Nptos=[];
for k=1:12
    mesk=k;
    day1=monthlims(k);
    day30=monthlims(k+1); %aprox. 30 days con dato para el mesk.
    Xk=X(:,day1:day30); %datos para el mesk.
    mXk=mean(Xk,2);
    stdXk=std(Xk,0,2);
    mX=[mX,mXk];
    stdX=[stdX,stdXk];
end
%..........
% $$$ figure(200)
% $$$ for k=1:12
% $$$     mXk=mX(:,k);
% $$$     subplot(3,4,k)
% $$$     plot(mXk,[1:m]*(-1))
% $$$ end
