function [Cb, CG, CXX, V, Sigma] =get_covSURFd(X,Y,b,indEq)

% N is the number of equations

resid  = Y- X*b;
N = length(unique(indEq));
n      = size(resid,1)/N;% number of observations per equation
K      = size(X,2)/N;
resid  = reshape(resid, n,N);

Sigma  = resid'*resid/n;% cov of error

V      = kron(Sigma, ones(K,K));

CXX    = cell(N,N);

for k=1:N
    indk  = indEq == k;
    Xk    = X(indk,((1:K)+(k-1)*K));
    for l=1:k
    indl  = indEq == l;
    Xl    = X(indl,((1:K)+(l-1)*K));
    CXX{k,l} = Xk'*Xl/n;  
    end
end

for k=2:N
    for l=1:k-1
    CXX{l,k} = CXX{k,l}';  
    end
end

CXX   = cell2mat(CXX);
CG     = V.*CXX;

invCXXDiag = inv(X'*X/n);% only within equations and not bewteen equations

Cb     = invCXXDiag*CG*invCXXDiag;

end






