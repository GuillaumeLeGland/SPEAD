function [H,VarH]=index_SaW(A,base)
%INDEX_SAW Shannon-Wiener index, 
%  also incorrectly known as Shannon-Weaver,is one of
%  several diversity indices used to measure versity
%  in categorical data. It is simply the Information
%  entropy of the distribution, treating species as
%  symbols and their relative population sizes as the
%  probability.
%  The advantage of this index is that it takes into
%  account the number of species and the evenness of
%  the species. The index is increased either by
%  having additional unique species, or by having a
%  greater species evenness. 
%  The value of the Shannon-Wiener Index usually lies 
%  between 1.5 and 3.5 for ecological data and rarely
%  exceeds 4.0. According to Southwood and Henderson
%  (2000), it is an insensitive measure of the character
%  of the S:N relationship and is dominated by the
%  abundant species.
%                                 s  
%                                 
%                           H= - sum{ pi * log(pi) }
%                                 
%                                 i=1
%  where, 
%        pi= ni/N, The relative abundance of each species.      
%        ni= The number of individuals in species i; the 
%            abundance of species i. 
%         N=  The total number of all individuals. 
%         s=  The number of species. Also called species 
%             richness.
%
%
%  Input: 
%         A= matrix of number of individuals of the species
%           (row=species and columns= stations of sample).
%           If no specie in the sample, put 0 or NaN.
%
%         base= base of the logarithm (default=2)
%
%  Output:
%         H =    index of Shannon-Weiner
%      VarH =    variance of the index of Shannon-Weiner
%
%
%  Example 1. From the Table 6.14 of number of individuals belonging to
%  different insect species obtained in collections using 20 light traps at 
%  the Nelliampathy location, taken from the FAO Corporate Document 
%  Repository: A Statistical Manual for Forestry Research
%  (http://www.fao.org/docrep/003/x6831e/X6831E14.htm)
%         A=[91,67,33,22,27,23,12,14,11,10,9,9,5,1,4,2,2,1,2,4];
%         [H,VarH]=index_SaW(A,exp(1))
%
%  Answer is:
%         H = 2.3717
%         VarH = 0.0029
%
%  Example 2. From the Table 6.14 of number of individuals belonging to
%  different insect species obtained in collections using 20 light traps at 
%  the Nelliampathy location, taken from the FAO Corporate Document 
%  Repository: A Statistical Manual for Forestry Research
%  (http://www.fao.org/docrep/003/x6831e/X6831E14.htm)
%
%    A=[91,84;
%       67,60;
%       33,40;
%       22,26;
%       27,24;
%       23,20;
%       12,16;
%       14,13;
%       11,12;
%       10,7;
%        9,5;
%        9,5;
%        5,9;
%        1,4;
%        4,6;
%        2,2;
%        2,4;
%        1,4;
%        2,5;
%        4,1];
%
%  Answer is:
%            H = 2.3717    2.4484
%         VarH = 0.0029    0.0027
%
%   Created by R. Hernandez-Walls and A. Trujillo-Ortiz 
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Ensenada, Baja California
%             Mexico.
%             rwalls@uabc.edu.mx
%
%   To cite this file, this would be an appropriate format:
%   Hernandez-Walls, R.  and Trujillo-Ortiz, A. (2010).
%   index_SaW:Shannon-Wiener Index. A MATLAB file.
%   http://www.mathworks.com/matlabcentral/fileexchange/28591
%
%   References:
%   Southwood, T. R. E. and Henderson, P. A. (2000), Ecological Methods 3rd
%      Edition. Blackwell Science 575pp 
%   Wikipedia, Shannon index. (2010), 
%      http://en.wikipedia.org/wiki/Shannon_index
%

if nargin < 2
      base = 2; %default log2
end
[r,c] = size(A);
warning('off','MATLAB:dispatcher:InexactMatch')
if r==1
    A=A';
end
A(find(isnan(A))) = 0;
n = sum(A); 
N = repmat(n,r,1);
B = A./N.*log(A./N)./log(base);
B(find(isnan(B))) = 0;
H = -sum(B);

%variance of the Shannon-Wiener index
aa = A./A;
aa(find(isnan(aa))) = 0;
S = sum(aa);
VarH = (sum(N./A.*(B.^2))-(sum(B)).^2)./n+(S-1)./(2*n.^2);
