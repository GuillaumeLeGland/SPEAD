function  df = mynumjac( values, pt, fcn, varargin)
% numjac calcul numérique du jacobien de de la fonction fcn
%
% df = numjac( values, pt, fcn, d)
%
% paramètres d'entrée
% -> values 
%        valeur de fcn au point "pt"
% -> pt
%        point d'évaluation
% -> fcn
%        fonction d'un vecteur "a" de même dimension que "pt" et des 
%        arguments supplémentaires d
%        fcn doit être de la forme  vals = fcn(a, varargin)
%        et vals est un vecteur colonne.
%
% paramètres de sortie
% <- df
%        Jacobien de fcn au point pt
%
%<http://www.mathworks.com/matlabcentral/fileexchange/18404>

df = zeros(length(values), length(pt));
for j = 1:length(pt)
   temp =  pt(j);
   h =  sqrt(eps)*abs(temp); 
   if (h == 0.0)
      h = sqrt(eps);
   end;
   pt(j) = temp+h;
   h =  pt(j)-temp;
% $$$    f=feval(fcn,pt,varargin{:}); %Original.
   f=feval(fcn,0,pt); %for ode45.
   for i = 1:length(values)
      df(:,j)= (f-values)/h; 
   end;
   pt(j) = temp;
end;
