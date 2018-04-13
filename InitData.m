function [data, boptions,testMode] = InitData(fileName,boptions,testMode,dlength)

if ~testMode.flag
    varIndex = [1 2 7];
    oIndex = [6];
    
    X = [];y = [];
    iter = 1;
    readrow = sprintf('C%d:I%d',iter+3,iter+3);
    while 1
        Excel_data = xlsread(fileName,1,readrow);
        if(isempty(Excel_data) || iter==dlength+1)
            break;
        else
            X(iter,:) = Excel_data(varIndex);
            y(iter) = Excel_data(oIndex);
        end
        iter = iter + 1;
        readrow = sprintf('C%d:I%d',iter+3,iter+3);
    end
    
    load(boptions.SearchGridFileName);
    data.max_x = max(SearchGrid,[],1);
    multp = (pinv(diag(data.max_x)));
    SearchGrid = SearchGrid(1:225,:);
    SearchGrid = SearchGrid*multp;
    X = X*pinv(diag(data.max_x));
    data.X = X;
    [SearchGrid] = removeXfromGrid(SearchGrid,data.X);
    boptions.SearchGrid = SearchGrid;
    mean_y = 0;
    max_y = 10;
    min_y = 1;
    var_y = ((max_y-min_y)/sqrt(12))^2;
    data.mean_y = mean_y;
    data.var_y = var_y;
    y = (y-mean_y)/sqrt(data.var_y);
    data.y = y';
    
else
    
    nDim = testMode.nDim;
    bounds = testMode.bounds;
    if strcmp(fileName,'Source')
        [X,y] = sourceData(testMode);        
    elseif strcmp(fileName,'Target')
        N0 = nDim + 1;
        X = repmat(bounds(:,1)',N0,1)+repmat([bounds(:,2)-bounds(:,1)]',N0,1).*rand(N0,nDim);        
        [X,y] = targetData(testMode,X);
    else
        fprintf('Error - please specify function type\n');
    end
    data.min_x = min(bounds,[],2);
    data.max_x = max(bounds,[],2);
    boptions.bounds = bounds;
    data.X = X;
    mean_y = 0;
    max_y = 1;
    min_y = 0;
    var_y = 1;
    data.mean_y = mean_y;
    data.var_y = var_y;
    y = (y-mean_y)/sqrt(data.var_y);
    data.y = y;
end






