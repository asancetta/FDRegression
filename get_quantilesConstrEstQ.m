function q = get_quantilesConstrEstQ(CG,H, A,c,lb,ub, nG, p)

%nG is number of simulations used to derive quantile at pr level p

[V,D]         = eig(CG);
thresh        = 0;
D(D < thresh)   = 0;
CGRoot   = sqrt(real(D))*V';% square root of CG;

K        = size(CG,1);
G        = randn(nG,K)*CGRoot; % simulate the Gaussian vector 
dBetas   = zeros(nG,K);
Aeq    = [];
ceq    = [];

for l = 1:nG
    fl = G(l,:)';
    [dbl,fval,exitflag,output]          = quadprog(H,fl,A,c,Aeq,ceq,lb,ub);
    dBetas(l,:)          = dbl;
end

q    = quantile(dBetas, p)';
end