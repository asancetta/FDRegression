function [A, c] = get_constr(indVar,k,a)

nA                                = length(a)-1;
K1                                = length(indVar);
ind1                              = find(indVar == k);
k1                                = length(ind1);
A                                 = zeros(k1-nA,K1);
for i = 1:(k1-nA)
    
    A(i,ind1(i:i+nA))             = a;

end
c                                 = zeros(k1-nA,1);


end
