function [X, y] = targetData(testMode,X)

if(strcmp(testMode.mode,'synthetic'))
    nDim = testMode.nDim;
    muVec = testMode.targetMean(1)*ones(1,nDim);
    SigmaMat = (1^2)*eye(nDim);
    cost = mvnpdf(X,muVec,SigmaMat);
    y = 1 - cost*(sqrt((2*pi)^nDim*det(SigmaMat)));
elseif(strcmp(testMode.mode,'real'))
    
    xxtrain = cell2mat(testMode.xxtrain(testMode.targetDigit));
    yytrain = cell2mat(testMode.yytrain(testMode.targetDigit));
    xxtest = cell2mat(testMode.xxtest(testMode.targetDigit));
    yytest = cell2mat(testMode.yytest(testMode.targetDigit));
    
    nSamples = size(X,1);
    y = zeros(nSamples,1);
    
    for ii = 1:nSamples
        gamma = (10^X(ii,1));
        cost = 2^X(ii,2);
        libsvm_options = sprintf('-s 0 -t 2 -g %f -c %f -q',gamma, cost);
        model = svmtrain(yytrain, xxtrain, libsvm_options);
        [~, accuracy, ~] = svmpredict(yytest, xxtest, model, '-q');
        y(ii) = accuracy(1)/100;
    end
    y = 1 - y;
elseif(strcmp(testMode.mode,'elastic_net'))
    nSamples = size(X,1);
    y = zeros(nSamples,1);
    xxtrain = cell2mat(testMode.xxtrain(testMode.targetDigit));
    yytrain = cell2mat(testMode.yytrain(testMode.targetDigit));
    xxtest = cell2mat(testMode.xxtest(testMode.targetDigit));
    yytest = cell2mat(testMode.yytest(testMode.targetDigit));

    glmOpt.label = yytrain;        
    for ii = 1:size(X,1)
        %%glmBin options
        lambda1 = 10^X(ii,1);
        lambda2 = 10^X(ii,2);
        glmOpt.l2_penalty		= lambda2; % vary it from 1e-2 to 1e-5
        glmOpt.l1_penalty		= lambda1; %keep this as zero
        glmOpt.nIters			= 10000;
        glmOpt.epsilon			= 1e-9;
        glmOpt.report_interval	= 100;
        
        fitall = glmBin('logit',[xxtrain ones(size(xxtrain,1),1)],zeros(1,size(xxtrain,2)+1),glmOpt,'train'); %bias is also added
        %%test
        deci = glmBin('logit',[xxtest ones(size(xxtest,1),1)],fitall,glmOpt,'test');
        [~,~,~,auc] = perfcurve(yytest,deci,1);
        y(ii) = auc;
    end
    y = 1 - y;
else
    fprintf('Unsupported mode\n');  
end