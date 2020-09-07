function []=mymonthlyaxis()
    
%********************************************
%Programa MYMONTHLYAXIS.m: Este programa me de los ejes como J,F,M,etc.
%
%Use: mymonthlyaxis
%********************************************
monthlim=[1,31,59,90,120,151,181,212,243,273,304,334,365];
xticklabel=['J','F','M','A','M','J','J','A','S','O','N','D','J']';
set(gca,'Xlim',[1 365],'Xtick',monthlim,'XTickLabel',xticklabel,'FontSize',[4])
grid on
