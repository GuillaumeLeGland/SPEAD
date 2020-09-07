close all
clear all

A = rand(3,5);
isub = [1,3] %rows 
jsub = [1,3,5] %cols 
Alowres = ones(length(isub),length(jsub))*nan; 
[JSUB,ISUB] = meshgrid(jsub,isub); %(x-cols,y-rows)
[INDEX] = sub2ind(size(A),ISUB,JSUB);
Alowres = A(INDEX)



