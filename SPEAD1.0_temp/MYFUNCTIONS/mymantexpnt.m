function [arg_coef,arg_expo] = mymantexpnt(arg)
% MANTEXPNT
% Returns the mantissa and exponent of a real base 10 argument.
% FILENAME: mantexpnt.m
% 20070710

% TEST DATA:
% arg = (0.5 - rand) * 10^fix(10*(0.5-rand))
% sprintf('\n\t%23.15E\n',arg)

for j = 1:length(arg)
    argj = arg(j);
    sgn = sign(argj);
    expnt = fix(log10(abs(argj)));
    mant = sgn * 10^(log10(abs(argj))-expnt);
    if abs(mant) < 1
	mant = mant * 10;
	expnt = expnt - 1;
    end
    arg_expo(j) = expnt;
    arg_coef(j) = mant;
end

% sprintf('\n\t%6.4f x E%+04d\n',mant,expnt)
%*******************************************
return % Added to indicate end of function

arg = 0.032;
arg_exp = floor(log10(arg));
arg_coef = arg / 1*10^arg_exp;

