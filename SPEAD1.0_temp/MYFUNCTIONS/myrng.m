function [] = myrng(sd)
% <https://www.mathworks.com/help/matlab/math/updating-your-random-number-generator-syntax.html>
rand('twister',5489) %rng('default')
rand('state',sd) %rng(sd)
rand('seed',sd)  %rng(sd)
return
