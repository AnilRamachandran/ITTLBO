function [dataTargetBrd] = bardenet(dataSource,dataTargetBrd,testMode,paramTarget,boptions)

    [Xst_brd,yst_brd] = mergeSourceTarget(dataSource.X,dataSource.y,dataTargetBrd.Xt,dataTargetBrd.ytn);

    dataTargetBrd.X = Xst_brd;
    dataTargetBrd.y = yst_brd;
    dataTargetBrd.max_x = max(testMode.bounds,[],2);%1*ones(testMode.nDim,1);
    dataTargetBrd.min_x = min(testMode.bounds,[],2);%0*ones(testMode.nDim,1);
    dataTargetBrd.max_y = 1;
    dataTargetBrd.min_y = 0;
    dataTargetBrd.mean_y = 0;
    dataTargetBrd.var_y = 1;
    
    gpTarget2_brd = buildGP(dataTargetBrd,paramTarget,boptions);

% recommend new sample 
    xnew2_brd = recommendSampleTransfer(gpTarget2_brd);
  
    newdata = newdataPoint(gpTarget2_brd,testMode,xnew2_brd);

    dataTargetBrd.Xt = [dataTargetBrd.Xt;newdata.X];
    dataTargetBrd.y0 = [dataTargetBrd.y0;newdata.y]; 
    dataTargetBrd.ytn = rankedfunction(dataTargetBrd.Xt,dataTargetBrd.y0');
