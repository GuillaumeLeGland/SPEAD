close all
clear all

kw=0.04; %(kw) irradiance attenutation due to water [m-1]
kp=0.03; %(kp) irradiance attenutation due to phyto self-sheding [m2*mmolN-1]

P=ones(1,10)*0.2;

f = inline('-(kw+kp*P)');
q = quad(f,0,5)
