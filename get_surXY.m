function [X1,Y1, indVar, indEq] = get_surXY(X,Y,ind)
% group X and Y by index 

[nN, K]     = size(X);

indUnique  =  unique(ind);
    
indUnique  = get_columVec(indUnique)';
nIndUnique = length(indUnique);
n          = nN/nIndUnique;

X1         = zeros(n,nIndUnique*K);
Y1         = zeros(n,1);

indVar     = zeros(nIndUnique*K,1);
indEq      = zeros(nN,1);
s          = 1;
for l = indUnique

  blockl   = ind == l;  
  Y1(1+(n*(s-1)):n*s) = Y(blockl,:);
  X1(1+(n*(s-1)):n*s,1+K*(s-1):K*s)  = X(blockl,:);
  indVar(1+K*(s-1):K*s)     = (1:K)';
  indEq(1+(n*(s-1)):n*s)  = s;
  s   = s+1;
end



end