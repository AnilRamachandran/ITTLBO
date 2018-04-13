function xnew = recommendSampleTransfer(gp)

boptions = gp.boptions;
optMethod = boptions.optMethod;
kvec = 1.0*ones(1,size(gp.X,1));
xinit = zeros(1,size(gp.X,2));
OptTime = 1;
if strcmp(optMethod,'SearchGrid')
    SearchGrid = boptions.SearchGrid;
    gVal = zeros(1,size(SearchGrid,1));
    for gg = 1 : size(SearchGrid,1)
        gVal(gg) = acqFuncTransfer(SearchGrid(gg,:)',gp,boptions);
    end
    xnew = SearchGrid(gVal==min(gVal),:);
    pp = randperm(size(xnew,1));
    try
        xnew = xnew(pp(1),:);
    catch
        keyboard
    end
else
    ymin = min(gp.y);
    [maxf, xnew] = recommendSample(gp.X,gp.y,ymin,gp.N,kvec,gp.invK,...
                boptions.bounds(:,1),boptions.bounds(:,2),gp.M,OptTime,gp.kernelTypeNum,gp.param.kernelScale,...
                gp.param.kernelVar,gp.param.rqalpha,boptions.eps,gp.msrSigma2scalar,boptions.criteriaNum,...
                xinit,1);

end