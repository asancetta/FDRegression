function out = estimate_volume_ved()
%Except where stated otherwise, all functions are Copyright (c) Alessio
%Sancetta, can be freely reproduced and used, but I assumes no responsibility.
% This code produces the empirical results for the paper 
%"Intraday End of Day Volume Prediction" by Alessio Sancetta.
%This code has been tested and runs on Matlab 2014a and requires some toolboxes, 
%but dependency is minimal.
 
outDir           = '~/VED/Empirical/';%put output dir here'
outDirGraph      = [ outDir 'Graphs/'];
inDirData        = '~/VED/DATA/';%put data dir here




lineWidth   =1.6;
ccyPair     = { 'FUT_6E';'FUT_6J';'FUT_6S';'FUT_6B';'FUT_6A';'FUT_6N';'FUT_6C';'FUT_ES'};




if ~exist(outDir, 'dir')

    mkdir(outDir);

end

if ~exist(outDirGraph, 'dir')

    mkdir(outDirGraph);

end

             

%% download data
load([inDirData,'dataVolumes.mat']);

%%the data loaded are cumulative volumes at one minute frequency for each
%%instrument (ftMat), the epoch time in nanoseconds with day of the week (ref_time),  
%%and a vector of zeros with entries equal to one if the row corresponds to 
%%the end of day (endDayFlag) 

weekDay                           = ref_time(:,2);
ref_time                          = ref_time(:,1);

indVolumes  = (1:size(ftMat,2));

indInfraDay     = (1:length(unique(ref_time)))'; 
indInfraDay     = repmat(indInfraDay,sum(endDayFlag),1);


indDay     = cumsum(endDayFlag);
indDay     = get_lag(indDay,1);
indDay(1)  = 0;
indDay     = indDay+1;


%the variables to use
K             = length(indVolumes);
assert(length(ccyPair)== K)
volume        = cell(1,K);
dVolume       = cell(1,K);
volume2goL    = cell(1,K);
volumeL       = cell(1,K);
volume2go     = cell(1,K); 
volumeEnd     = cell(1,K);
volumeL5      = cell(1,K);
volumeLAvg    = cell(1,K);
volume2goLAvg = cell(1,K);
volume2goL5   = cell(1,K); 
X             = cell(1,K);
Y             = cell(1,K);
for k = 1:K

[volume{k}, dVolume{k}, volume2goL{k}, volumeL{k}, volume2go{k}, volumeEnd{k},...
                    volumeL5{k}, volumeLAvg{k},volume2goLAvg{k},volume2goL5{k}] = ...
                                  get_volumeVariables(ftMat, indDay,...
                                                       indVolumes(k), weekDay);

X{k}                              = [ones(size(volume{k},1),1), volume{k}, ...
                                     volumeL{k},volume2goL{k}+volumeL{k}];
                                 
Y{k}                              = volume2go{k};
end


%% plot a few sample realizations for the data
for k =1:K
    
close all
volumeFData  = get_fData(volume{k}, max(indInfraDay),0)';
pause(3)
plot(volumeFData(:,1:10));
set(findobj(gca,'Type','text'),'FontSize',12)
set(gca,'FontSize', 14)
ylabel(['Cumulative Volume']);
xlabel(['Time in Minutes']);
set(findobj(gca, 'Type', 'Line'), 'LineWidth', lineWidth);

graphName           = 'volumePlot';
outFile           =[outDirGraph  graphName,ccyPair{k}, '.png']; 
saveas(gcf, outFile);  

end
%compute summary stats and save
summaryStats   =  cell(K,1);
summaryStats2   =  cell(K,1);

for k =1:K
    volumeFData  = get_fData(volume{k}, max(indInfraDay),0)';
    dataMask     = diff(volumeFData,1);
    mF           = mean(dataMask,2);
    dataMask2    = bsxfun(@rdivide, dataMask, mF);
    dataMask     = vec(dataMask);
    resAutocorr  = autocorr(dataMask,120);
    [resSummary, colNames] = get_summaryStats(dataMask);
    summaryStats{k}     = [resSummary([1:4]),resAutocorr([6,101])', resSummary(end-1)];
    %--------repeat the same but with deseasonalised data 
    dataMask     = vec(dataMask2);
    resAutocorr  = autocorr(dataMask,120);
    [resSummary, colNames] = get_summaryStats(dataMask);
    summaryStats2{k}     = [resSummary([1:4]),resAutocorr([6,101])', resSummary(end-1)];


end

summaryStats  = cell2mat(summaryStats);
summaryStats  = num2cell(summaryStats);
cols          = [{'mean'}, {'std'},{'skew'}, {'kurt'},{'acf5'},{'acf100'},{'n'}] ;
out           = [cols;summaryStats];

summaryStats2  = cell2mat(summaryStats2);
summaryStats2  = num2cell(summaryStats2);
cols          = [{'mean_deas'}, {'std_deas'},{'skew_deas'}, {'kurt_deas'},{'acf5_deas'},{'acf100_deas'},{'n'}] ;
out           = [out;cols;summaryStats2];

outFile           =  'summaryStats.csv'; 
writecell2csv(out ,outDir, outFile);



%% delete data at start and end of day and nans:
indDel    = (indInfraDay == min(indInfraDay)) | (indInfraDay == max(indInfraDay));
for k =1:K
    indDel = indDel |  any(isnan([Y{k},X{k}]),2);
    
end

        
for k = 1:K
    X{k}(indDel,:) = [];
    Y{k}(indDel,:) = [];
    volumeEnd{k}(indDel,:) = [];
end

indInfraDay(indDel,:) = [];
indDay(indDel,:)      = [];
XES                   = X{K};
%% loop over ccy's
resMSE1 = [];
resMSE2 = [];

PRIdaily =[];

for k =1:(K-1)
        suffix                 = ccyPair{k};  
       [resMSETmp,PRIdailyTmp] =  get_estimatedByCcy(X{k},Y{k},indDay,indInfraDay,XES,volumeEnd{k},...
                                                suffix, outDir, outDirGraph);
    
       resMSE1 = [resMSE1;resMSETmp([2],:)];
       resMSE2 = [resMSE2;resMSETmp([3],:)];
    
       PRIdaily =[PRIdaily;PRIdailyTmp];
    
end
resMSEHeading = resMSETmp([1],:);
resMSE =[resMSE1;resMSE2]; 
csvwrite([outDir, 'PRI_Daily_All','.csv'],PRIdaily) 

writecell2csv([resMSEHeading;resMSE] ,outDir, ['PRI_MSE_All','.csv']);

out     = 1;
toc
end



%% Estimation functions

%estimation wrapper for each instrument
function [resMSE,PRIdaily] = get_estimatedByCcy(X,Y,indDay,indInfraDay,X6A, volumeEnd, suffix, outDir, outDirGraph)

%%
%indModelPresentation = [6,1,2,3,4,5];
indIn            = ismember(indDay, (1:70));
indOut           = ~indIn;

lineWidth        = 1.6;

if size(X6A,1)>0;
    X6Ain   = X6A(indIn,:);
else
    X6Ain   = [];
end

simpleFlag               = 0; %uses the one step constrained estimator
[b1,b2, b3,b4, b6A,b16A] = get_estimatorVolumes(X(indIn,:),Y(indIn,:),...
                                                indInfraDay(indIn,:),X6Ain,simpleFlag);

[pred1, pred2, pred3, pred4, pred5,pred6] = get_predVolumes(X,  indInfraDay,...
                                                          b1,b2,b3,b4, b6A,b16A, X6A);

%% residuals diagnostics
error   = bsxfun(@minus, Y, [pred1,pred2,pred3,pred4, pred5, pred6]);
labelsPlot = {'H1','H2','H3','H4','H5','U'};
nScores    = 3; 
nLags      = 5;
autocorrError(error, indInfraDay, nScores, nLags, labelsPlot,outDir, outDirGraph, suffix );
%% PRI and MSE
mY2     = mean(Y.^2);
mean(abs(error(indIn,:)).^2)
std(abs(error(indIn,:)).^2)/sqrt(sum(indIn))

mean(abs(error(indOut,:)).^2)
std(abs(error(indOut,:)).^2)/sqrt(sum(indOut))


MSE  = mean(abs(error(indOut,:)).^2)
PRI  = 100*(1-bsxfun(@rdivide, MSE(1:end-1),MSE(end)))
PRI1A  = 100*(1-bsxfun(@rdivide, MSE(2:end),MSE(1)));

PRI1C  = 100*(1-bsxfun(@rdivide, MSE(3),MSE(5)));


% test unconstrained versus fully constrained and unconstrained versus
% constrained augmented
errorOut2 = error(indOut,:).^2;
indDayOut = indDay(indOut);
uniqueDays  = unique(indDayOut);
nDayOut   = size(uniqueDays,1);
test      = zeros(nDayOut,size(errorOut2,2));
for s=1:nDayOut
    
    indMask   = indDayOut == uniqueDays(s);
    test(s,:)  =   mean(errorOut2(indMask,:))/mY2;
    
end

testA     = bsxfun(@minus, test(:,2:end), test(:,1));
testB    = bsxfun(@minus, test(:,1:end-1), test(:,end));

test1A     = sum(testA)./sqrt(sum(testA.^2))
test1B     = sum(testB)./sqrt(sum(testB.^2))

testC    = bsxfun(@minus, test(:,3), test(:,5))
test1C     = sum(testC)./sqrt(sum(testC.^2))


%% graph

close all
pause(3)
labelsBoxPlot = labelsPlot;
boxplot(test,'labels',labelsBoxPlot, 'datalim', [0,3*max(iqr(test))] );
set(findobj(gca,'Type','text'),'FontSize',12)
set(gca,'FontSize', 14)
ylabel(['Square Error']);
set(findobj(gca, 'Type', 'Line'), 'LineWidth', lineWidth);

graphName           = 'boxPlotSquareError';
outFile           =[outDirGraph  graphName suffix '.png']; 
saveas(gcf, outFile);  

%% compare to simple estimator just for illustration
simpleFlag                     = 1; %uses the two step constrained estimator: not that good
[b1s,b2s, b3s,b4s, b6As,b16As] = get_estimatorVolumes(X(indIn,:),Y(indIn,:),...
                                                      indInfraDay(indIn,:),X6Ain,simpleFlag);

                                
                                
[pred1s, pred2s, pred3s, pred4s, pred5s, pred6s] = get_predVolumes(X, indInfraDay,...
                                                          b1s,b2s,b3s,b4s, b6As,b16As, X6A);

errors   = bsxfun(@minus, Y, [pred1s,pred2s,pred3s,pred4s, pred5s,pred6s]);

mean(abs(errors(indIn,:)).^2)
std(abs(errors(indIn,:)).^2)/sqrt(sum(indIn))

mean(abs(errors(indOut,:)).^2)
std(abs(errors(indOut,:)).^2)/sqrt(sum(indOut))

MSEs  = mean(abs(errors(indOut,:)).^2)
PRIs  = 100*(1-bsxfun(@rdivide, MSEs(1:end-1),MSEs(end)))

PRI1As  = 100*(1-bsxfun(@rdivide, MSEs(2:end),MSEs(1)));

PRI1Cs  = 100*(1-bsxfun(@rdivide, MSEs(3),MSEs(5)));

errorOut2 = errors(indOut,:).^2;
indDayOut = indDay(indOut);
uniqueDays  = unique(indDayOut);
nDayOut   = size(uniqueDays,1);
testS      = zeros(nDayOut,size(errorOut2,2));
for s=1:nDayOut
    
    indMask    = indDayOut == uniqueDays(s);
    testS(s,:) =   mean(errorOut2(indMask,:))/mY2;
    
end

testSA     = bsxfun(@minus, testS(:,2:end), testS(:,1));
testSB    = bsxfun(@minus, testS(:,1:end-1), testS(:,end));

test2A     = sum(testSA)./sqrt(sum(testSA.^2))
test2B     = sum(testSB)./sqrt(sum(testSB.^2))


testSC    = bsxfun(@minus, testS(:,3), testS(:,5));
test2C     = sum(testSC)./sqrt(sum(testSC.^2))

resMSE     = [PRI ,PRIs; test1B, test2B];


%% augment results and add cols and row names
resMSE     = [resMSE,[PRI1A,PRI1As; test1A, test2A]];

resMSE     = [resMSE,[PRI1C,PRI1Cs; test1C, test2C]];


colsResMSE = [{''}, strcat(labelsPlot(1:5),['vsU']),strcat(labelsPlot(1:5),['vsU_S'])];
colsResMSE = [colsResMSE, strcat(labelsPlot(2:6),['vsH1']),strcat(labelsPlot(2:6),['vsH1_S'])];
colsResMSE = [colsResMSE, strcat(labelsPlot(3),['vsH5']),strcat(labelsPlot(3),['vsH5_S'])];

rowsResMSE = strcat(suffix,{'_PRI', '_tStat'})';

resMSE     = num2cell(resMSE);
resMSE     = [colsResMSE;[rowsResMSE,resMSE]];

%% daily predictions
maxLag   = 8;

Vdaily   = volumeEnd(indInfraDay == max(indInfraDay));
intercept = ones(size(Vdaily));
Xdaily   = [intercept, lagmatrix(Vdaily,1:maxLag)]; 
Xdaily(isnan(Xdaily)) = 0;
for l = 1:maxLag
    Xtmp         = Xdaily(:,1:1+l);
    beta         = (Xtmp'*Xtmp)\(Xtmp'*Vdaily);
    res{l}       = beta;
    
end

%use the largest model
predDaily = kron(Xdaily*res{end}, ones(size(unique(indInfraDay))) );
predIntraDay = pred1+X(:,2);

MSEdaily   = mean([(volumeEnd-predDaily),(max(volumeEnd-predDaily,0)),(volumeEnd-pred1-X(:,2))].^2);
PRIdaily1 = 100*(1-bsxfun(@minus, MSEdaily(2:end),MSEdaily(1)));
PRIdaily2 = 100*(1-bsxfun(@minus, MSEdaily(end),MSEdaily(2)));

PRIdaily = [PRIdaily1,PRIdaily2]; 


%% put together all results and save

writecell2csv(resMSE ,outDir, ['PRI_MSE',suffix,'.csv']);
csvwrite([outDir, 'PRI_Daily', suffix,'.csv'],PRIdaily) 

%% generate graphs 

indZoom  = ismember(indDay, [123:127]); 

data2plot =  [volumeEnd,predDaily,predIntraDay];
colors  = {'r','g','b'};
legendNames = {'Actual','Daily','IntraDay'};
close all

for l=1:size(data2plot,2)
    plot(data2plot(indZoom,l),'color',colors{l})
    hold on
end

legend(legendNames, 'location', 'NorthWest')
ylabel('Volumes');
xlabel('Time in Minutes');  
set(findobj(gca,'Type','text'),'FontSize',12)
set(gca,'FontSize', 14)
set(findobj(gca, 'Type', 'Line'), 'LineWidth', lineWidth);


graphName           = ['volumePredictionLast5Days'];
outFile           =[outDirGraph  graphName suffix '.png']; 
saveas(gcf, outFile);  


%% transform data to construct confidence bands
bu   = b3;
[XS,YS, indVar, indEq] = get_surXY(X(indIn,:),Y(indIn,:),indInfraDay(indIn,:));
[Cb, CG, CXX] = get_covSURFd(XS,YS,bu, indEq);
varBu = diag(Cb);

indDayIn = indDay(indIn);
uniqueDaysIn  = unique(indDayIn);
nDaysIn   = size(uniqueDaysIn,1);


indk = [1,2,3,4]; 
for k= indk
    close all
    pause(3)
    indk = (indVar==k); 
    sdBu = sqrt(varBu(indVar ==k))/sqrt(nDaysIn);
    qu   = [-1.96*sdBu,1.96*sdBu];
    bandsu  = bsxfun(@plus, bu(indk), qu);
    plot(b1(indk),'color','b');
    hold on
    plot(b3(indk),'color','r');
    hold on
    plot(bandsu(:,1), 'linestyle', '--');
    hold on
    plot(bandsu(:,2), 'linestyle', '--');
    set(findobj(gca,'Type','text'),'FontSize',12)
    set(gca,'FontSize', 14)
    ylabel(['b',sprintf('%d', k)]);
    xlabel('Time in Minutes');  
    set(findobj(gca, 'Type', 'Line'), 'LineWidth', lineWidth);

    graphName           = ['b_' sprintf('%d', k)];
    outFile           =[outDirGraph  graphName  suffix '.png']; 
    saveas(gcf, outFile);  
end    
 
%% get confidence bands for estimator under binding constraint for b4
%tangent cone: convex decreasing (4)
a      = [-1,1];
k      = 4;
[A1, c1] = get_constr(indVar,k,a);

a      = -[1,-2,1];
k      = 4;
[A2, c2] = get_constr(indVar,k,a);
A2  =[];c2=[];
A      = [A1;A2];
c      = [c1;c2];
H      = XS'*XS;
ub     = 200*ones(size(XS,2),1);
lb     = -ub;
qc     =  get_quantilesConstrEstQ(CG,H/nDaysIn, A,c,lb,ub, 500, [.025, .975]);
qc1    = [qc(indVar==k,1),qc(indVar==k,2)]/sqrt(nDaysIn);
bands1c  = bsxfun(@plus, mean(bu(indVar==k,1)), qc1);
close all
plot(b1(indVar==k,1), 'color', 'b')
hold on;
plot(bu(indVar==k,1),'color', 'r')
hold on
plot(bandsu(:,1), 'linestyle', '--');
hold on
plot(bandsu(:,2), 'linestyle', '--');
hold on
plot(bands1c(:,1), 'linestyle', '-.','color', 'r');
hold on
plot(bands1c(:,2), 'linestyle', '-.','color', 'r');
set(findobj(gca,'Type','text'),'FontSize',12)
set(gca,'FontSize', 14)
ylabel(['b',sprintf('%d', k)]);
xlabel('Time in Minutes');
set(findobj(gca, 'Type', 'Line'), 'LineWidth', lineWidth);
graphName           = ['b_b0Boundary' sprintf('%d', k)];
outFile           =[outDirGraph  graphName suffix '.png']; 
saveas(gcf, outFile);  

close all

end

%computes and saves residuals diagnostics
function autocorrError(error, indInfraDay, nScores, nLags, labelsPlot, outDir, outDirGraph, suffix )


[n,K]= size(error);
Lmin = min(indInfraDay);
Lmax = max(indInfraDay);
L    = Lmax-Lmin+1;

outDirGraphACF  = [outDirGraph, 'ACFS/'];

if ~exist(outDirGraphACF, 'dir')

    mkdir(outDirGraphACF);

end


chis_df   = nScores*nLags;
statsACF  = ones(K,2);
for k = 1:K
    errorMat         = reshape(error(:,k),L,n/L)';
    covError         =  cov(errorMat);
    [errorV, errorD] = eig(covError);
    [eigenvalues, indSort]  = sort(diag(errorD),'descend');
    errorV                  = errorV(:,indSort); 

    
    
    errorS           =  errorMat*errorV(:,1:nScores);
    resACF           = zeros(nLags, nScores);
    for s = 1:nScores
        autocorr(errorS(:,s))
        set(findobj(gca,'Type','text'),'FontSize',14)
        set(gca,'FontSize', 14)
        ylabel('ACF');
        xlabel('Lag');
        title([labelsPlot{s} ': Factor Score ', sprintf('%d', s)] )
        graphName           = ['ACF_Resid_Score',sprintf('%d', s),'_',labelsPlot{k},'_'];
        outFile           =[outDirGraphACF, graphName  suffix '.png']; 
        saveas(gcf, outFile);  
        close all
        
        acfTmp       = autocorr(errorS(:,s),nLags);
        resACF(:,s)  = acfTmp(2:nLags+1);           
    end
    statsACF(k,1)      = (n/L)*sum(vec(resACF.^2));
    statsACF(k,2)      = 1-chis_cdf(statsACF(k,1),chis_df);
end


summaryStats  = num2cell(statsACF);
cols          = {'BoxPierceStat','p-value'} ;
out           = [cols;summaryStats];
out           = [['Hypothesis';labelsPlot'],out];
outFile       =  ['acfStats', suffix, '.csv']; 
writecell2csv(out ,outDir, outFile);




end

function F = chis_cdf (x, a)
% PURPOSE: returns the cdf at x of the chisquared(n) distribution
%---------------------------------------------------
% USAGE: cdf = chis_cdf(x,n)
% where: x = a vector
%        n = a scalar parameter
% NOTE: chis_cdf(x,n) = gamm_cdf(x/2,n/2)
%---------------------------------------------------
% RETURNS:
%        a vector pdf at each element of x from chisq(n) distribution      
% --------------------------------------------------
% SEE ALSO: chis_d, chis_pdf, chis_rnd, chis_inv
%---------------------------------------------------

%        Anders Holtsberg, 18-11-93
%        Copyright (c) Anders Holtsberg
% documentation modified by LeSage to
% match the format of the econometrics toolbox

if (nargin ~= 2)
    error ('Wrong # of arguments to chis_cdf');
end

if any(any(a<=0))
   error('chis_cdf: dof is wrong')
end

F = gamm_cdf(x/2,a*0.5);
end

function cdf = gamm_cdf (x, a)
% PURPOSE: returns the cdf at x of the gamma(a) distribution
%---------------------------------------------------
% USAGE: cdf = gamm_cdf(x,a)
% where: x = a vector 
%        a = a scalar gamma(a)
%---------------------------------------------------
% RETURNS:
%        a vector of cdf at each element of x of the gamma(a) distribution      
% --------------------------------------------------
% SEE ALSO: gamm_d, gamm_pdf, gamm_rnd, gamm_inv
%---------------------------------------------------

%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg

if nargin ~= 2
error('Wrong # of arguments to gamm_cdf');
end;

if any(any(a<=0))
   error('gamm_cdf: parameter a is wrong')
end

cdf = gammainc(x,a);
I0 = find(x<0);
cdf(I0) = zeros(size(I0));
end



%% carries out the constrained estimations
function [b1,b2, b3, b4, b6A,b16A, indVar,indVar6A] = get_estimatorVolumes(X,Y, indInfraDay, X6A, simpleFlag)



[X0, indVar] = get_surX(X,indInfraDay);

[n,K]                             = size(X);

H            = X0'*X0;
f            = -X0'*Y;

%% unconstrained estimator

Aeq  = [];
ceq  = [];
A    = [];
c    = [];
lb   = [];
ub   = [];

[b1,fval,exitflag,output]          = quadprog(H,f,A,c,Aeq,ceq,lb,ub);


if simpleFlag>0
    H            = eye(size(b1,1));
    f            = -b1;
end

%% constrained 1:  decreasing  intercept

%decreasing intercept
a      = [-1,1];
k      = 1;
[A1, c1] = get_constr(indVar,k,a);

Aeq  = [];
ceq  = [];
A    = A1;
c    = c1;
lb   = -200*ones(1,size(X0,2));
ub   = 200*ones(1,size(X0,2));


[b2,fval,exitflag,output]          = quadprog(H,f,A,c,Aeq,ceq,lb,ub);

%% constrained 2: convex and concave decreasing (1), convex decreasing (2),...
%% convex decreasing(3), decreasing (4)


%convex and concave decreasing (1)
a      = [-1,1];
k      = 1;
[A1, c1] = get_constr(indVar,k,a);

indVar1 = zeros(size(indVar)); 
indVar2 =  zeros(size(indVar));
midPoint = round(find(indVar == 1, 1,'last')/2); 
indVar1(1:midPoint) = indVar(1:midPoint); 
indVar2(1+midPoint:end) = indVar(1+midPoint:end); 
a1     = -[1,-2,1];%convex
a2     = [1,-2,1];%concave
k1     = 1;k2 = 1;
[A2, c2] = get_constrMixed(indVar1,k1,a1,indVar2,k2,a2);

%convex decreasing (2)

a      = [-1,1];
k      = 2;
[A3, c3] = get_constr(indVar,k,a);

a      = -[1,-2,1];
k      = 2;
[A4, c4] = get_constr(indVar,k,a);


% %convex decreasing (3)

a      = [-1,1];
k      = 3;
[A5, c5] = get_constr(indVar,k,a);

a      = -[1,-2,1];
k      = 3;
[A6, c6] = get_constr(indVar,k,a);


% decreasing (4)
a      = [-1,1];
k      = 4;
[A7, c7] = get_constr(indVar,k,a);

Aeq  = [];
ceq  = [];
A    = [A1;A2;A3;A4;A5;A6;A7];
c    = [c1;c2;c3;c4;c5;c6;c7];
lb   = -200*ones(1,size(X0,2));
ub   = 200*ones(1,size(X0,2));

[b3,fval,exitflag,output]          = quadprog(H,f,A,c,Aeq,ceq,lb,ub);

%% constrained 3:  convex concave decreasing (1), constant coefficient (2,3,4)


%convex concave decreasing (1)
a      = [-1,1];
k      = 1;
[A1, c1] = get_constr(indVar,k,a);

indVar1 = zeros(size(indVar)); 
indVar2 =  zeros(size(indVar));
midPoint = round(length(indVar)/2);
indVar1(1:midPoint) = indVar(1:midPoint); 
indVar2(1+midPoint:end) = indVar(1+midPoint:end); 
a1     = -[1,-2,1];%convex
a2     = [1,-2,1];%concave
k1     = 1;k2 = 1;
[A2, c2] = get_constrMixed(indVar1,k1,a1,indVar2,k2,a2);


% constant coefficient (2,3,4)

a      = [-1,1];
kSet   = [2,3,4];
[A3, c3] = get_constrSet(indVar,kSet,a);

Aeq  = A3;
ceq  = c3;
A    = [A1;A2];
c    = [c1;c2];
lb   = -200*ones(1,size(X0,2));
ub   = 200*ones(1,size(X0,2));

[b4,fval,exitflag,output]          = quadprog(H,f,A,c,Aeq,ceq,lb,ub);

if size(X6A,1) >0

%% augmented model


% without intercept 
[X06A, indVar6A] = get_surX(X6A(:,2),indInfraDay);

indVar6A        = [indVar; max(indVar)+indVar6A];
X1              = [X0,X06A];

[n,K]                             = size(X);

H            = X1'*X1;
f            = -X1'*Y;

b16A         = H\(-f);

if simpleFlag>0
    H            = eye(size(b16A,1));
    f            = -b16A;
end


% % decreasing (1,2,3,4,5,6,7)
a      = [-1,1];
kSet   = [1,2,3,4,5];%,6,7];
[A1, c1] = get_constrSet(indVar6A,kSet,a);

% convex (2,3)
a      = -[1,-2,1];
kSet   = [2,3];%,6,7];
[A2, c2 ] = get_constrSet(indVar6A,kSet,a);


% convex and concave (1)
endPoint  = find(indVar == 1, 1,'last');

indVar1 = zeros(size(indVar6A)); 
indVar2 =  zeros(size(indVar6A));
midPoint = round(length(indVar)/2);
indVar1(1:midPoint) = indVar(1:midPoint); 
indVar2(1+midPoint:endPoint) = indVar(1+midPoint:endPoint); 
a1     = -[1,-2,1];%convex
a2     = [1,-2,1];%concave
k1     = 1;k2 = 1;
[A3, c3] = get_constrMixed(indVar1,k1,a1,indVar2,k2,a2);



Aeq  = [];
ceq  = [];
A    = [A1;A2;A3];
c    = [c1;c2;c3];

lb   = -200*ones(1,size(X1,2));
ub   = 200*ones(1,size(X1,2));

[b6A,fval,exitflag,output]          = quadprog(H,f,A,c,Aeq,ceq,lb,ub);

else
    
    b6A  = [];
    b16A = [];
    indVar6A = [];
end

end


%% computes predictions
function [pred1, pred2, pred3, pred4, pred5, pred6] = get_predVolumes(X, indInfraDay,...
                                                             b1,b2,b3,b4, b5,b6, X6A)


[X0, indVar] = get_surX(X,indInfraDay);

pred1 = X0*b1;
pred2 = X0*b2;
pred3 = X0*b3;
pred4 = X0*b4;

if size(X6A,1)>0
[X06A, indVar6A] = get_surX(X6A(:,2),indInfraDay);

%indVar6A        = [indVar; max(indVar)+indVar6A];
X1             = [X0,X06A];
pred5           = X1*b5;
pred6           = X1*b6;
else
    pred5           = ones(size(pred3));
    pred6           = pred5;
end    
    
end

%% simulates to compute quantiles to construct confidence bands for estimator under binding constraints 
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

%% utility functions to compute constraints
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

function [A, c] = get_constrMixed(indVar1,k1,a1,indVar2,k2,a2)


[A1, c1] = get_constr(indVar1,k1,a1);
[A2, c2] = get_constr(indVar2,k2,a2);

A  = [A1;A2];
c  = [c1;c2];

end

function [A, c] = get_constrSet(indVar,kSet,a)

A      = [];
c      = [];
kSet   = get_columVec(kSet)';

for k = kSet

    [AMask, cMask] = get_constr(indVar,k,a);

    A   = [A;AMask];
    c   = [c;cMask];

end

end

%% transform variables for SUR
function [X1, indVar] = get_surX(X,ind)


[n, K]     = size(X);

indUnique  =  unique(ind);
    
indUnique  = get_columVec(indUnique)';
nIndUnique = length(indUnique);

X1         = zeros(n,nIndUnique*K);
indVar     = zeros(nIndUnique*K,1);
s          = 1;
for l = indUnique

  blockl   = ind == l;  
  X1(blockl,1+K*(s-1):K*s)  = X(blockl,:);
  indVar(1+K*(s-1):K*s)     = (1:K)';
  s   = s+1;
end



end

function [Cb, Csigma, CXX] =get_covSURFd(X,Y,b,indEq)

% N is the number of equations

resid  = Y- X*b;
N = length(unique(indEq));
n      = size(resid,1)/N;% number of observations per equation
K      = size(X,2)/N;

CXX    = cell(N,N);
Csigma = cell(N,N);


for k=1:N
    indk  = indEq == k;
    Xk    = X(indk,((1:K)+(k-1)*K));
    XkError   = bsxfun(@times, X(indk,((1:K)+(k-1)*K)), resid(indk));
    
    for l=1:k
    indl  = indEq == l;
    Xl    = X(indl,((1:K)+(l-1)*K));
    XlError   = bsxfun(@times, X(indl,((1:K)+(l-1)*K)), resid(indl));
    CXX{k,l} = Xk'*Xl/n;  
    Csigma{k,l} = XkError'*XlError/n;  
    
    end
end


for k=2:N
    for l=1:k-1
    CXX{l,k} = CXX{k,l}';
    Csigma{l,k} = Csigma{k,l}';
    end
end


CXX    = cell2mat(CXX);
Csigma = cell2mat(Csigma);

invCXXDiag = inv(X'*X/n);% only within equations and not bewteen equations

Cb     = invCXXDiag*Csigma*invCXXDiag;

end

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


function X = get_columVec(X)
% it also works with a matrix, in which case, it transposes if the nCol is
% larger than nRow

if size(X,1)< size(X,2)
   
    X        = X';
    
end



end

%% data parsers utils

%produce the variables to be used
function [volume, dVolume, volume2goL, volumeL, volume2go, volumeEnd,...
          volumeL5, volumeLAvg, volume2goLAvg, volume2goL5] = ...
                                  get_volumeVariables(ftMat1, indDay,indVol, weekDays)
                              
blocks          = unique(indDay)';
volume          = ftMat1(:,indVol);
volume2goL      = zeros(size(ftMat1,1),1);
dVolume         = zeros(size(ftMat1,1),1);
volumeL         = zeros(size(ftMat1,1),1);
volumeLTmp      = [];
volumeEnd       = zeros(size(ftMat1,1),1);
volume2go       = zeros(size(ftMat1,1),1);

for i = blocks
   
    indTmp      = indDay == i;
    volumeTmp   = volume(indTmp,:);
    if isempty(volumeLTmp)
       volumeLTmp   = volumeTmp;
    end
    
    volumeEnd(indTmp,:) = max(volumeTmp);
    volume2go(indTmp,:) = max(volumeTmp)-volumeTmp;
    volume2goL(indTmp,:)= max(volumeLTmp)-volumeLTmp;
    volumeL(indTmp,:)   = volumeLTmp;
    dVolume(indTmp,:)   = volumeTmp - volumeLTmp;
    volumeLTmp          = volumeTmp;
    
end


N = sum(indDay==blocks(1));
nDays  = max(indDay);

uniqueWeekDays = unique(weekDays)';

volumeFd    = reshape(volume,N,nDays )';
volumeL5    = zeros(nDays,N);

volume2goFd    = reshape(volume2go,N,nDays )';
volume2goL5    = zeros(nDays,N);

weekDaysFd  = reshape(weekDays,N,nDays )';

assert(all(all(diff(weekDaysFd') ==0)), 'Sampling problem with FData')
weekDaysFd    = weekDaysFd(:,1);
for s = uniqueWeekDays
    
    indWeek   = weekDaysFd == s;
    volumeL5Mask  = volumeFd(indWeek,:);
    volumeL5(indWeek,:) = get_lag(volumeL5Mask,1);

    volume2goL5Mask  = volume2goFd(indWeek,:);
    volume2goL5(indWeek,:) = get_lag(volume2goL5Mask,1);

end
volumeL5   = vec(volumeL5');

volume2goL5   = vec(volume2goL5');



volumeLFd   = reshape(volumeL,N,nDays )';
volumeLAvg   = nanmean(volumeLFd,2)./volumeLFd(:,end);
volumeLAvg   = bsxfun(@times, volumeLAvg, ones(nDays,N));
volumeLAvg   = vec(volumeLAvg');

volume2goLFd    = reshape(volume2goL,N,nDays )';
volume2goLAvg   = nanmean(volume2goLFd,2)./volumeLFd(:,end);
volume2goLAvg   = bsxfun(@times, volume2goLAvg, ones(nDays,N));
volume2goLAvg   = vec(volume2goLAvg');



end

%reshape data to functional data format
function [out, index]= get_fData(Y,K, startFlag)

if nargin <3
    
   startFlag  =1;
    
end


if startFlag ==1
    N            = floor((size(Y,1)-1)/K);
    index        = zeros(N,K);
    jStart       = 2;
    jEnd         = K+1;
    startVal     = Y(1,:);
    out          = zeros(N,K);

    for i =1:N
            index(i,:)       = jStart:jEnd;
            out(i,:)       = Y(jStart:jEnd,:)-startVal;
            startVal       = Y(jEnd,:);
            jStart         = jStart+K;
            jEnd           = jEnd+K;
    end

else
    N            = floor(size(Y,1)/K);
    index        = zeros(N,K);
    jStart       = 1;
    jEnd         = K;
    out          = zeros(N,K);

    for i =1:N
            index(i,:)       = jStart:jEnd;
            out(i,:)       = Y(jStart:jEnd,:);
            jStart         = jStart+K;
            jEnd           = jEnd+K;
    end


    
end

    
end

% produce summary stats
function [out, colNames] = get_summaryStats(x)

indNan       = isnan(x);

x(indNan,:)  = [];

out(1)       = mean(x);
out(2)       = std(x);
out(3)       = skewness(x);
out(4)       = kurtosis(x);
out(5)       = min(x);
out(6)       = max(x);
out(7)       = median(x);
out(8)       = iqr(x);
out(9)       = length(indNan);
out(10)      = sum(indNan);

colNames     = {'mean', 'std', 'skew' , 'kurt', 'min','max', ...
                'med', 'iqr', 'n', 'NaN'};

end

function output = get_lag(data, lag)
%% computes the lagged matrix of data or leading matrix by shift and inserting NaN's; use negative lag for lead 
K = size(data,2);


if lag>=0
output = [repmat(NaN, lag,K); data(1:size(data,1)-lag, :)];
else
lag = -lag ;    
output = [   data( (1+lag):size(data,1),:  ); repmat(NaN,lag,K)    ];
end

output;
end

function [y]=vec(x)
% PURPOSE:
% Creates a column vector by appending the columns of a matrix to each other. 
% 
% FORMAT:
% yc = vec( x); 
%
% INPUT:
% x     NxK matrix. 
% 
% OUTPUT 
% yc   (N*K)x1 vector, the columns of x appended to each other. 
% 
% REMARKS:
% 
% EXAMPLE: 
% 
%      x = [ 1 2,
%            3 4 ];
%      yc = vec(x);
% 
% x =    1.000000    2.000000 
%        3.000000    4.000000 
%   
% yc =    1.000000 
%         3.000000 
%         2.000000 
%         4.000000 
%
% Written by Dimitris Korobilis (2007-2008)
% University of Strathclyde
  
y=reshape(x,size(x,1)*size(x,2),1);

end


%% write to csv's utils:
function writecell2csv(out, outDir, fileName)


if ~exist(outDir, 'dir')
    mkdir(outDir);
end

outFile           = [outDir, fileName];    
    
fid = fopen(outFile , 'w');

%fprintf(fid , '%s\n' , ''); 
writeCSV(fid, out);
fprintf( fid, '\n');

fclose(fid);    
end


function writeCSV(fid, data)


[N, K]  = size(data);

spc = cell(1,K);
for k=1:K-1
   spc{k}   = ',';  
end
spc{K}      = '';

for n = 1:N
    for k=1:K
        if isnumeric(data{n,k}) 
            dataString   = sprintf('%f', data{n,k});
        else
            dataString   = data{n,k};
        end
        dataString   = [dataString, spc{k}];
        fprintf(fid , '%s' , dataString);

    end
    fprintf( fid, '\n');
end

end

%%
