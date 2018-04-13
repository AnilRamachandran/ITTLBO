function [X, y] = sourceData(testMode)

fprintf('Building source\n');
if(strcmp(testMode.mode,'synthetic'))
    nDim = testMode.nDim;
    N0 = 50;
    bounds = testMode.bounds;
    X = repmat(bounds(:,1)',N0,1)+repmat([bounds(:,2)-bounds(:,1)]',N0,1).*rand(N0,nDim);
    muVec = testMode.sourceMean(1)*ones(1,nDim);
    SigmaMat = (1^2)*eye(nDim);
    cost = 2*(mvnpdf(X,muVec,SigmaMat));
    y = (1 - cost*(sqrt((2*pi)^nDim*det(SigmaMat))));
elseif(strcmp(testMode.mode,'real'))
    bounds = testMode.bounds;
    g1 = [bounds(1,1):bounds(1,2)];
    g2 = [bounds(2,1):bounds(2,2)];
    [xg1,xg2] = meshgrid(g1,g2);
    X = [xg1(:) xg2(:)];
    nSamples = size(X,1);
    y = zeros(nSamples,1);
    
    xxtrain = cell2mat(testMode.xxtrain(testMode.sourceDigit));
    yytrain = cell2mat(testMode.yytrain(testMode.sourceDigit));
    xxtest = cell2mat(testMode.xxtest(testMode.sourceDigit));
    yytest = cell2mat(testMode.yytest(testMode.sourceDigit));
     
    for ii = 1:nSamples
        fprintf('%2d/%2d',ii,nSamples);
        gamma = (10^X(ii,1));
        cost = 2^X(ii,2);
        libsvm_options = sprintf('-s 0 -t 2 -g %f -c %f -q',gamma, cost);
        model = svmtrain(yytrain, xxtrain, libsvm_options);
        [~, accuracy, ~] = svmpredict(yytest, xxtest, model, '-q');
        y(ii) = accuracy(1)/100;
        fprintf('\b\b\b\b\b');
    end
    y = 1 - y;
elseif(strcmp(testMode.mode,'elastic_net'))
    addpath('glmBin');
    bounds = testMode.bounds;
    g1 = [bounds(1,1):bounds(1,2)];
    g2 = [bounds(2,1):bounds(2,2)];
    [xg1,xg2] = meshgrid(g1,g2);
    X = [xg1(:) xg2(:)];
    nSamples = size(X,1);
    y = zeros(nSamples,1);
    
    xxtrain = cell2mat(testMode.xxtrain(testMode.sourceDigit));
    yytrain = cell2mat(testMode.yytrain(testMode.sourceDigit));
    xxtest = cell2mat(testMode.xxtest(testMode.sourceDigit));
    yytest = cell2mat(testMode.yytest(testMode.sourceDigit));
    
    for ii = 1:nSamples
        fprintf('%2d/%2d',ii,nSamples);
        
        %%glmBin options
        lambda1 = 10^X(ii,1);
        lambda2 = 10^X(ii,2);
        glmOpt.l2_penalty		= lambda2; % vary it from 1e-2 to 1e-5
        glmOpt.l1_penalty		= lambda1; %keep this as zero
        glmOpt.nIters			= 10000;
        glmOpt.epsilon			= 1e-9;
        glmOpt.report_interval	= 100;
        
        glmOpt.label = yytrain;
        fitall = glmBin('logit',[xxtrain ones(size(xxtrain,1),1)],zeros(1,size(xxtrain,2)+1),glmOpt,'train'); %bias is also added
        %%test
        deci = glmBin('logit',[xxtest ones(size(xxtest,1),1)],fitall,glmOpt,'test');
        [~,~,~,auc] = perfcurve(yytest,deci,1);
        y(ii) = auc;
        fprintf('\b\b\b\b\b');
    end
    y = 1 - y;
else
    fprintf('Unsupported mode\n');
end