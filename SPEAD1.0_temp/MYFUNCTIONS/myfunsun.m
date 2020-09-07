function [y]=myfunsun(ymax,ymin,par,parmax,parmin)

%**********************************************************************
%Program MYFUNSUN.m: Este programa genera una variable "y=f(par)" de
%size(365) que varia entre 'ymin' e 'ymax' (eq. de LEFEVRE etal. 2002).
%
%Use: [y]=myfunsun(ymax,ymin,par,parmax,parmin)
%
%par: vector(365) or array(depth,365) con la irradiacion solar.
%**********************************************************************
%-------------------------------------------------------------------
%NOTA: Si ymin=0 y parmin=0 esta ecuacion se convierte en: ymax*(par/parmax) (eg. GABRIC'05, TellusB 57)
%-------------------------------------------------------------------

y = ymax - (ymax-ymin)*((parmax-par)./(parmax-parmin));
