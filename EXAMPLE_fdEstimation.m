% EXAMPLE OF WORK FLOW

N   = 50; % number of samples in a day
n   = 100; % number of days

%% auxiliary variables
time = (1:N)/N;% fraction of time of day 
indInfraDay = repmat((1:N)',n,1); % will tell where in the day we are in
indDay      = kron((1:n)',ones(N,1)); % will tell which day we are in

%% simulate some data
rng('default');
% an intercept
X1   = ones(n,N);


% smooth random trigonometric curve

X2   = zeros(n,N);
for j=1:3
    rng(j);
   
    X2   = X2+(2^(1-j))*bsxfun(@times,-log(rand(n,1)), cos(2*pi*j*time+(j/2)));
    
end
X2    = X2/6;

error   = zeros(n,N);
for j=1:3
    rng(100+2*j);

    error   = error+(2^(1-j))*bsxfun(@times,-log(rand(n,1)), cos(2*pi*j*time+(j/2)));
    
end
error  = .5*error/6;% a smooth error


b01     = (time.^3)'; % a convex increasing function;
b02     = ones(N,1);

X      = [vec(X1'),vec(X2')];% instructive to plot X1',X2' to see what they are 
Xb     = bsxfun(@times, X1,b01')+bsxfun(@times, X2,b02');
Y      = Xb+error;
% turn Y into a vector 
Y      = vec(Y');

% X Y are a vector and a matrix where each eantry is a different time of a
% day for different days. Use indDay to group by day and indInfraDay to group by time of day 

%% construct variables for estimation and perform estimation

% construct the variables for SUR regression, but uses quad programming
[XS,YS, indVar, indEq] = get_surXY(X,Y,indInfraDay);
%[XS, indVar] = get_surX(X,indInfraDay); %indInfraDay indentifies an equation
% the output indVar identifies the position in XS of each variable, 
% e.g. XS(:,(indVar == 1)) extract the first variable

H            = XS'*XS;
f            = -XS'*YS;

% peforms unconstrained estimation using quadratic programming

Aeq  = [];
ceq  = [];
A    = [];
c    = [];
lb   = [];
ub   = [];

[bu,fval,exitflag,output]          = quadprog(H,f,A,c,Aeq,ceq,lb,ub);

%extract bu for each of two variables
bu1   = bu(indVar == 1); bu2   = bu(indVar == 2); 

%plots results for each variable separately
close all
figure(1); 
pause(1);
subplot(2,1,1); plot([b01,bu1]);
title(['b1 - Estimated and True'] )
subplot(2,1,2); plot([b02,bu2]);
title(['b2 - Estimated and True'] )

% perform constrained estimation:
% constraint is convex increasing intercept (variable 1)

k   = 1; % index of variable 1
a   = [1,-1]; % constraint for increasing variable
[A1 , c1] = get_constr(indVar,k,a);

a   = -[1,-2,1]; % constraint for convex variable
[A2 , c2] = get_constr(indVar,k,a);

% constraint is A*b<=c 
A   = [A1;A2]; 
c   = [c1;c2];

[bc,fval,exitflag,output]          = quadprog(H,f,A,c,Aeq,ceq,lb,ub);

%extract bc for each of the two variables
bc1   = bc(indVar == 1); bc2   = bc(indVar == 2); 

% plot results 

figure(2); 
pause(1);
subplot(2,1,1); plot([b01,bc1]);
title(['b1 - Estimated and True'] )
subplot(2,1,2); plot([b02,bc2]);
title(['b2 - Estimated and True'] )

%% compare error with true function

db2  = [b01-bu1,b02-bu2,b01-bc1,b02-bc2].^2;
mean(db2)
std(db2)

% relative improvement (in MSE) when using a constraint 

1-mean(db2(:,3:4))./mean(db2(:,1:2))% for b1 and b2

%% get confidence bands for estimators

[Cb, CG, CXX, V] =get_covSURFd(XS,YS,bu, indEq);

varBu = diag(Cb);

sdBu1 = sqrt(varBu(indVar ==1))/sqrt(n);
sdBu2 = sqrt(varBu(indVar ==2))/sqrt(n);

qu1   = [-1.96*sdBu1,1.96*sdBu1];
qu2   = [-1.96*sdBu2,1.96*sdBu2];

bands1u  = bsxfun(@plus, bu1, qu1);
bands2u  = bsxfun(@plus, bu2, qu2);

%Observe that the first 17 values of bc1 are constant. 
% You decide to draw confidence bands assuming that the true value
% is constant over the first 17 values. This is equivalent to say that the constraint
%is binding over the first 17 points: 
% convex increasing constraint is binding only if the function is constant.
% The tanget cone at the true (hypothesized) value is 
% such that difference of the first 17 coefficients is less than zero
% and the second difference of the first 17 coefficients is less than zero:

A   = [A1(1:16,:);A2(1:15,:)];
c   = [c1(1:16);c2(1:15)];

qc     =  get_quantilesConstrEstQ(CG,H/n, A,c,lb,ub, 500, [.025, .975]);

qc1    = [qc(indVar==1,1),qc(indVar==1,2)]/sqrt(n);
qc2    = [qc(indVar==2,1),qc(indVar==2,2)]/sqrt(n);

bands1c  = bsxfun(@plus, bc1, qc1);
bands2c  = bsxfun(@plus, bc2, qc2);


figure(3); 
pause(1);
subplot(2,1,1); plot([b01,bu1,bands1u]);
title(['b1 - Estimated with Confidence Bands'] )
subplot(2,1,2); plot([b02,bu2,bands2u]);
title(['b2 - Estimated with Confidence Bands'] )


figure(4); 
pause(1);
subplot(2,1,1); plot([b01,bc1,bands1c]);
title(['b1 - Estimated with Confidence Bands'] )
subplot(2,1,2); plot([b02,bc2,bands2c]);
title(['b2 - Estimated with Confidence Bands'] )
% observe smoother confidence bands for bc1 under the hypothesis that
% constraint is binding for true b01 function for first 17 points. 
%%  perform constrained estimation directly from the unconstrained estimator

Hs   = eye(length(bu));
fs   = -bu;

k   = 1; % index of variable 1
a   = [1,-1]; % constraint for increasing variable
[A1 , c1] = get_constr(indVar,k,a);

a   = -[1,-2,1]; % constraint for convex variable
[A2 , c2] = get_constr(indVar,k,a);

% constraint is A*b<=c 
A   = [A1;A2]; 
c   = [c1;c2];

[bcs,fval,exitflag,output]          = quadprog(Hs,fs,A,c,Aeq,ceq,lb,ub);

%extract bc for each of three variables
bcs1   = bcs(indVar == 1); bcs2   = bcs(indVar == 2); 

% plot results 

figure(5); 
pause(1);
subplot(2,1,1); plot([b01,bcs1]);
title(['b1 - Estimated and True'] )
subplot(2,1,2); plot([b02,bcs2]);
title(['b2 - Estimated and True'] )

%% compare error with true function

db2  = [b01-bu1,b02-bu2,b01-bcs1,b02-bcs2].^2;
mean(db2)

% relative improvement (in MSE) when using constraint 
1-mean(db2(:,3:4))./mean(db2(:,1:2))% for b1 and b2
% this is not as good as the estimator bc.
