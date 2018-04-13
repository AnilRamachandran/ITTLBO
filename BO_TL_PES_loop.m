% clear;
addpath ../GPstuff-4.4/diag
addpath ../GPstuff-4.4/dist
addpath ../GPstuff-4.4/gp
addpath ../GPstuff-4.4/mc
addpath ../GPstuff-4.4/misc
addpath ../GPstuff-4.4/optim
addpath ../GPstuff-4.4/xunit
addpath ../sourceFiles

displayFlag = 1;
testMode.flag = 1;

if testMode.flag==1
    testMode.nDim = 3;
    fileNameSource = 'Source';
    fileNameTarget = 'Target';
    testMode.delta = 10*0.97;
    testMode.mode = 'elastic_net';
    if(strcmp(testMode.mode,'synthetic'))
        testMode.targetMean = [0.3];
        srcMeanVec = [-1.8;-0.7;0.4;1.5];
        bounds = repmat([-2 2],testMode.nDim,1);
    elseif(strcmp(testMode.mode,'real'))
        testMode.nDim = 2;
        bounds(1,:) = [-3 3];
        bounds(2,:) = [-3 3];
        testMode = read_classfi_dataasets(testMode);
        testMode.sourceDigits = [2,3,4,5];
        load svm_target_digit
        testMode.targetDigit = svm_target_digit;
     elseif(strcmp(testMode.mode,'elastic_net'))
        testMode.nDim = 2;
        bounds = repmat([-2 0],testMode.nDim,1);
        testMode = read_classfi_dataasets(testMode);
        testMode.sourceDigits = [2,3,4,5];
        load en_target_digit
        testMode.targetDigit = en_target_digit;
    else        
    end
    testMode.bounds = bounds;
    
else
    fprintf('\nerror in specifying run mode');
end
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
boptions.eps = 0.0;
boptions.optMethod = 'Continuous';%'SearchGrid', 'Continuous'

%% Set Parameters
kernelType = 'SE';%'RQ', 'SE'
if strcmp(kernelType,'SE')
    kernelTypeNum = 1;
elseif strcmp(kernelType,'RQ')
    kernelTypeNum = 2;
else
    kernelTypeNum = 1;
end

paramSource.kernelType = kernelType;
paramSource.kernelScale = 0.1;%*sqrt(testMode.nDim);
paramSource.kernelVar = 1;
paramSource.rqalpha = 10;
paramSource.msrSigma2 = (0.01)^2;
paramSource.kernelTypeNum = kernelTypeNum;

paramTarget.kernelType = kernelType;
paramTarget.kernelScale = 0.1;%*sqrt(testMode.nDim);
paramTarget.kernelVar = 1;
paramTarget.rqalpha = 10;
paramTarget.msrSigma2 = (0.01)^2; % measurement noise variance
paramTarget.kernelTypeNum = kernelTypeNum;


if(strcmp(testMode.mode,'synthetic'))
    [dataTarget, boptions] = InitData(fileNameTarget,boptions,testMode);
    NumSrc = size(srcMeanVec,1);
    for ss = 1 : NumSrc
        testMode.sourceMean = srcMeanVec(ss,:);
        [dataSource(ss), boptions] = InitData(fileNameSource,boptions,testMode);
        gpSource(ss) = buildGP(dataSource(ss),paramSource,boptions);
    end
elseif(strcmp(testMode.mode,'real'))
    NumSrc = length(testMode.sourceDigits);
    for ss = 1 : NumSrc
        testMode.sourceDigit = testMode.sourceDigits(ss);
        [dataSource(ss), boptions] = InitData(fileNameSource,boptions,testMode); 
        gpSource(ss) = buildGP(dataSource(ss),paramSource,boptions);
    end
    [dataTarget, boptions] = InitData(fileNameTarget,boptions,testMode);
elseif(strcmp(testMode.mode,'elastic_net'))
    NumSrc = length(testMode.sourceDigits);
    for ss = 1 : NumSrc
        testMode.sourceDigit = testMode.sourceDigits(ss);
        [dataSource(ss), boptions] = InitData(fileNameSource,boptions,testMode); 
        gpSource(ss) = buildGP(dataSource(ss),paramSource,boptions);  
    end
    [dataTarget, boptions] = InitData(fileNameTarget,boptions,testMode);
end
dataSource_merged.y = [];dataSource_merged.X = [];
for ii = 1:length(dataSource)
    dataSource_merged.X = [dataSource_merged.X;dataSource(ii).X];
    dataSource_merged.y = [dataSource_merged.y;dataSource(ii).y];
    dataSource_merged.min_x = dataSource(ii).min_x;
    dataSource_merged.max_x = dataSource(ii).max_x;
    dataSource_merged.mean_y = dataSource(ii).mean_y;
    dataSource_merged.var_y = dataSource(ii).var_y;
end
gpSourceMerged = buildGP(dataSource_merged,paramSource,boptions);

dataSource_mergedbrd.X = [];dataSource_mergedbrd.y = [];
for ii = 1:length(dataSource)
    yt = [];
    dataSource_mergedbrd.X = [dataSource_mergedbrd.X;dataSource(ii).X];
    yt = rankedfunction(dataSource(ii).X,dataSource(ii).y');
    dataSource_mergedbrd.y = [dataSource_mergedbrd.y;yt'];
    dataSource_mergedbrd.min_x = dataSource(ii).min_x;
    dataSource_mergedbrd.max_x = dataSource(ii).max_x;
    dataSource_mergedbrd.mean_y = dataSource(ii).mean_y;
    dataSource_mergedbrd.var_y = dataSource(ii).var_y;
end

%% Build Target GP (No transfer)
gpTarget = buildGP(dataTarget,paramTarget,boptions);

dataTargetT2 = dataTarget;
dataTarget_MergedSrc = dataTarget;
dataTargetYog = dataTarget;
dataTargetBrd = dataTarget;

dataTargetYog.Xt = dataTargetYog.X;
dataTargetYog.y0 = dataTargetYog.y;
dataTargetYog.y_dummy = dataTargetYog.y;


dataTargetBrd.Xt = dataTargetBrd.X;
dataTargetBrd.y0 = dataTargetBrd.y;
dataTargetBrd.ytn = rankedfunction(dataTargetBrd.X,dataTargetBrd.y');

%% Generating samples from the posterior
simulation = load('simulation.txt');
num_Samples_transfer = 100; %200 elastic net
Rndmfeatures = 500;%1000
num_Samples_notransfer = num_Samples_transfer;
xstar = {};
for ss = 1:NumSrc    
    % We sample from the global minimum using sources
    [ l_src, sigma_src, sigma0_src ] = sampleHypers(dataSource(ss).X, dataSource(ss).y, num_Samples_transfer);
	[ m_src hessians_src] = sampleMinimum(num_Samples_transfer, dataSource(ss).X, dataSource(ss).y, sigma0_src, sigma_src, l_src, testMode.bounds(:,1), testMode.bounds(:,2), Rndmfeatures);
    xstar{ss}.sample = m_src;
    xstar{ss}.hessian = hessians_src;
end


%% Bayesian optimization loop
iter_t = 1;% required for GP-UCB
opt_iter = 1;
Maxiters = 100;
guesses = dataTarget.X;
guessesT2 = dataTargetT2.X;
percentage_matrix = [];
[ l, sigma, sigma0 ] = sampleHypers(dataTargetT2.X, dataTargetT2.y, num_Samples_notransfer);
lT2 = l;
sigmaT2 = sigma;
sigma0T2 = sigma0;

lT2_mat = [];sigmaT2_mat = [];sigma0T2_mat = [];l_mat = [];sigma_mat = [];sigma0_mat = [];
for opt_iter = 1 : Maxiters
%% Transfer learning method with Source x_star samples
    
    % We sample from the global minimum using target
    xstar{NumSrc+1} = [];
	[ m_tar hessians_tar] = sampleMinimum(num_Samples_transfer, dataTargetT2.X, dataTargetT2.y, sigma0T2, sigmaT2, lT2, testMode.bounds(:,1), testMode.bounds(:,2), Rndmfeatures);
    xstar{NumSrc+1}.sample = m_tar;
    xstar{NumSrc+1}.hessian = hessians_tar;
    
%% Divergence between source and target global minimum distributions
    dist_mat = [];   
    distanceMeasure = 'KL'; 
    if strcmp(distanceMeasure,'KL')
         for kk = 1:NumSrc
            distNNsrc = [];
            for  ii=1:num_Samples_transfer
                distNNsrc_temp = [];
                for jj=1:num_Samples_transfer
                    if ii == jj
                        continue;
                    else
                        distNNsrc_temp = [distNNsrc_temp;norm((xstar{kk}.sample(ii)-xstar{kk}.sample(jj)),2)];
                    end
                end
                distNNsrc = [distNNsrc;min(distNNsrc_temp)];
            end
            distNNsrc2tar = [];
            for  ii=1:num_Samples_transfer
                distNNsrc2tar_temp = [];
                for jj=1:num_Samples_transfer
                    distNNsrc2tar_temp = [distNNsrc2tar_temp;norm((xstar{kk}.sample(ii)-xstar{NumSrc+1}.sample(jj)),2)];
                end
                distNNsrc2tar = [distNNsrc2tar;min(distNNsrc2tar_temp)];
            end
            dist = (testMode.nDim/num_Samples_transfer)*sum(log(distNNsrc2tar./distNNsrc))+log(num_Samples_transfer/(num_Samples_transfer-1));
            dist_mat = [dist_mat dist];
         end
    end
%%
    dist_mat = [dist_mat 0];
    similarity_measure = exp(-dist_mat/10);
    percentage_sample = (similarity_measure/sum(similarity_measure));
    percentage_matrix = [percentage_matrix percentage_sample'];
    
    mT2 = []; hessiansT2 = [];
    for i=1:length(percentage_sample)
        if percentage_sample(i) == 0
            continue;
        else
            numelements = round(percentage_sample(i)*length(xstar{i}.sample));
            indices = randperm(length(xstar{i}.sample));
            indices = indices(1:numelements); 
            mT2 = [mT2;xstar{i}.sample(indices,:)];
            hessiansT2 = [hessiansT2 xstar{i}.hessian(indices)];
        end
    end
    if size(mT2,1) < num_Samples_transfer
        indices = randperm(length(xstar{end}.sample));
        indices = indices(1:(num_Samples_transfer-size(mT2,1)));
        mT2 = [mT2;xstar{end}.sample(indices,:)];
        hessiansT2 = [hessiansT2 xstar{end}.hessian(indices)];
    else
        mT2 = mT2(1:num_Samples_transfer,:);
        hessiansT2 = hessiansT2(1:num_Samples_transfer);
    end

    
% We call the ep method
	retT2 = initializeEPapproximation(dataTargetT2.X, dataTargetT2.y, mT2, lT2, sigmaT2, sigma0T2, hessiansT2);
    
% We define the cost function to be optimized
	costT2 = @(x) evaluateEPobjective(retT2, x);
    
% We optimize globally the cost function
    optimumT2 = globalOptimizationOneArgument(costT2, testMode.bounds(:,1), testMode.bounds(:,2), guessesT2);

% Update TransferT2
    newdataT2 = newdataPoint(dataTargetT2,testMode,optimumT2);
    dataTargetT2 = augmentData(dataTargetT2, newdataT2);
    
% We sample from the posterior distribution of the hyper-parameters
    [ lT2, sigmaT2, sigma0T2 ] = sampleHypers(dataTargetT2.X, dataTargetT2.y, num_Samples_transfer);
     lT2_mat = [lT2_mat lT2];sigmaT2_mat = [sigmaT2_mat sigmaT2];sigma0T2_mat = [sigma0T2_mat sigma0T2];

% We update the kernel matrix on the samples
	KernelMatrixInvT2 = {};
	for j = 1 : num_Samples_transfer
		KernelMatrixT2 = computeKmm(dataTargetT2.X, lT2(j,:)', sigmaT2(j), sigma0T2(j));
        KernelMatrixInvT2{ j } = invChol_mex(KernelMatrixT2);
    end
    fT2 = @(x) posteriorMean(x, dataTargetT2.X, dataTargetT2.y, KernelMatrixInvT2, lT2, sigmaT2);
	gfT2 = @(x) gradientPosteriorMean(x, dataTargetT2.X, dataTargetT2.y, KernelMatrixInvT2, lT2, sigmaT2);
    
% We optimize the posterior mean of the GP		
	optimumT2 = globalOptimization(fT2, gfT2, testMode.bounds(:,1), testMode.bounds(:,2), guessesT2);
    guessesT2 = [ guessesT2 ; optimumT2 ];
    
% Display results
    if(testMode.flag)
        ymin = min(dataTargetT2.y);
        q1 = ymin*sqrt(dataTargetT2.var_y) + dataTargetT2.mean_y;
        q2 = newdataT2.y*sqrt(dataTargetT2.var_y) + dataTargetT2.mean_y;
        if displayFlag == 1
        fprintf('Transfer PES: ymin: %f, ycurrent: %f\n',q1,q2);
        end
    end
    
%% Transfer Learning Env-GP
% setting gp-ucb beta_t
    if strcmp(boptions.criteria,'gpUCB')
        delta_gpucb = 0.1;%0.01
        dd = testMode.nDim;
        aa = 0.1;%0.001;
        bb = 0.2;%1/2;
        z2 = pi*pi/6; %zeta(2)
        boptions.eps = 4*dd*log(dd*iter_t*bb*1*sqrt(log(2*dd*aa/delta_gpucb))) + ...
            2*log(4/delta_gpucb) + 2*(log(z2) + log(iter_t)*2);
        iter_t = iter_t + 1;
    end
    gpTarget.boptions.eps = boptions.eps;

    dataTargetENV.X = [dataSource_merged.X;dataTarget_MergedSrc.X];
    dataTargetENV.y = [dataSource_merged.y;dataTarget_MergedSrc.y];
    dataTargetENV.max_x = dataSource_merged.max_x;
    dataTargetENV.min_x = dataSource_merged.min_x;
    dataTargetENV.mean_y = dataSource_merged.mean_y;
    dataTargetENV.var_y = dataSource_merged.var_y;

    gpTargetENV = buildGP(dataTargetENV,paramTarget,boptions);
    gpTargetENV.alpha0 = 0.2;
    gpTargetENV.beta0 = 0.2;
    gpTargetENV.msrSigma2_source = (0.01)^2;
    gpTargetENV.Ns = size(gpSourceMerged.X,1);
    gpTargetENV.ys = gpSourceMerged.y;
    
    gpTargetENV = recomputeGPTargetT2(gpSource,gpTargetENV);% recomputing invK and yt on Xs
    xnewENV = recommendSampleTransfer(gpTargetENV);
    
    newdataENV = newdataPoint(gpTargetENV,testMode,xnewENV);
    gpTargetENV = addDataGP(gpTargetENV, newdataENV);
    dataTarget_MergedSrc.X = [dataTarget_MergedSrc.X;newdataENV.X];
    dataTarget_MergedSrc.y = [dataTarget_MergedSrc.y;newdataENV.y];  
   
    
% Display results
    if(testMode.flag)
        ymin = min(gpTargetENV.y(gpTargetENV.Ns+1:end));
        q1 = ymin*sqrt(gpTargetENV.var_y) + gpTargetENV.mean_y;
        q2 = newdataENV.y*sqrt(gpTargetENV.var_y) + gpTargetENV.mean_y;
        if displayFlag == 1
        fprintf('Transfer ENV-GP: ymin: %f, ycurrent: %f\n',q1,q2);
        end
    end
    
% Transfer Learning - Yogatama

    [dataTargetYog] = yogatama(dataSource,dataTargetYog,testMode,paramTarget,boptions);

% Display results
    if(testMode.flag)
        ymin = min(dataTargetYog.y0);
        q1 = ymin*sqrt(dataTargetYog.var_y) + dataTargetYog.mean_y;
        q2 = dataTargetYog.y0(end)*sqrt(dataTargetYog.var_y) + mean(dataTargetYog.mean_y);
        if displayFlag == 1
        fprintf('Transfer SMBO: ymin: %f, ycurrent: %f\n',q1,q2);
        end
    end

%% Transfer Learning - Bardnet
    [dataTargetBrd] = bardenet(dataSource_mergedbrd,dataTargetBrd,testMode,paramTarget,boptions); 
    
% Display results
    if(testMode.flag)
        ymin = min(dataTargetBrd.y0);
        q1 = ymin*sqrt(dataTargetBrd.var_y) + dataTargetBrd.mean_y;
        q2 = dataTargetBrd.y0(end)*sqrt(dataTargetBrd.var_y) + mean(dataTargetBrd.mean_y);
        if displayFlag == 1
        fprintf('Transfer SCoT: ymin: %f, ycurrent: %f\n',q1,q2);
        end
    end
    
%% No Transfer
% We sample from the global minimum using target
    [ m hessians] = sampleMinimum(num_Samples_notransfer, dataTarget.X, dataTarget.y, sigma0, sigma, l, testMode.bounds(:,1), testMode.bounds(:,2), Rndmfeatures); 
    
% We call the ep method
	ret = initializeEPapproximation(dataTarget.X, dataTarget.y, m, l, sigma, sigma0, hessians);
        
% We define the cost function to be optimized
	cost = @(x) evaluateEPobjective(ret, x);

% We optimize globally the cost function
    optimum = globalOptimizationOneArgument(cost, testMode.bounds(:,1), testMode.bounds(:,2), guesses);

% Update No transfer
    newdata = newdataPoint(dataTarget,testMode,optimum);
    dataTarget = augmentData(dataTarget, newdata);
    
% We sample from the posterior distribution of the hyper-parameters
    [ l, sigma, sigma0 ] = sampleHypers(dataTarget.X, dataTarget.y, num_Samples_notransfer);
    l_mat = [l_mat l];sigma_mat = [sigma_mat sigma];sigma0_mat = [sigma0_mat sigma0];
    
% We update the kernel matrix on the samples
	KernelMatrixInv = {};
	for j = 1 : num_Samples_notransfer
		KernelMatrix = computeKmm(dataTarget.X, l(j,:)', sigma(j), sigma0(j));
        KernelMatrixInv{ j } = invChol_mex(KernelMatrix);
    end
    f = @(x) posteriorMean(x, dataTarget.X, dataTarget.y, KernelMatrixInv, l, sigma);
	gf = @(x) gradientPosteriorMean(x, dataTarget.X, dataTarget.y, KernelMatrixInv, l, sigma);
    
% We optimize the posterior mean of the GP		
	optimum = globalOptimization(f, gf, testMode.bounds(:,1), testMode.bounds(:,2), guesses);
    guesses = [ guesses ; optimum ];
    
% Display results
    if(testMode.flag)
        ymin = min(dataTarget.y);
        q1 = ymin*sqrt(dataTarget.var_y) + dataTarget.mean_y;
        q2 = newdata.y*sqrt(dataTarget.var_y) + dataTarget.mean_y;
        if displayFlag == 1
        fprintf('Generic-BO(No transfer): ymin: %f, ycurrent: %f\n',q1,q2);   
        end
    end
end