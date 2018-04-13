function [yr] = rankedfunction(X,y)

%% Load Search Grid and set Bayesian Optimization Parameters

boptions.criteria = 'EI'; %'EI', 'PI','gpUCB'
if strcmp(boptions.criteria,'gpUCB')
    criteriaNum = 3;
elseif strcmp(boptions.criteria,'EI')
    criteriaNum = 1;
elseif strcmp(boptions.criteria,'PI')
    criteriaNum = 2;
else
    criteriaNum = 1;
end
boptions.criteriaNum = criteriaNum;

%% Set Parameters

kernelType = 'SE';%'RQ', 'SE'
if strcmp(kernelType,'SE')
    kernelTypeNum = 1;
elseif strcmp(kernelType,'RQ')
    kernelTypeNum = 2;
else
    kernelTypeNum = 1;
end

paramTarget.kernelType = kernelType;
paramTarget.kernelScale = 0.1;%*sqrt(testMode.nDim);
paramTarget.kernelVar = 1;
paramTarget.rqalpha = 10;
paramTarget.msrSigma2 = (0.01)^2; % measurement noise variance
paramTarget.kernelTypeNum = kernelTypeNum;


K = kernel(X,X,paramTarget);

if(length(paramTarget.msrSigma2)<2)
    gp.msrSigma2 = paramTarget.msrSigma2*ones(size(K,1),1);
else
    gp.msrSigma2 = paramTarget.msrSigma2;
end

gp.K = K+diag(gp.msrSigma2);

[gp.y, gp.invK] = fit_rankGP(y, gp.K, paramTarget);

gp.y1 = y;
yr = gp.y;


