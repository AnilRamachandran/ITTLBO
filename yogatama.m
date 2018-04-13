function [dataTargetYog] = yogatama(dataSource,dataTargetYog,testMode,paramTarget,boptions)

    dataSource_merged = [];
    for ii = 1:length(dataSource)
        [dataSource(ii).mean,dataSource(ii).var] = find_mean_var(dataSource(ii).y,'source');
        dataSource(ii).y_computed = compute_response(dataSource(ii).y,dataSource(ii).mean,dataSource(ii).var,'source');
        dataSource_new = [dataSource(ii).X dataSource(ii).y_computed];
        dataSource_merged = [dataSource_merged;dataSource_new];
    end
    
    [dataTargetYog.mean,dataTargetYog.var] = find_mean_var(dataTargetYog.y_dummy,'target');
    dataTargetYog.y_compute = compute_response(dataTargetYog.y_dummy,dataTargetYog.mean,dataTargetYog.var,'target');
    dataTargetYog_new = [dataTargetYog.Xt dataTargetYog.y_compute];
    
    dataMerged = [dataSource_new;dataTargetYog_new];

    dataTargetYog.X = dataMerged(:,1:end-1);
    dataTargetYog.Ns = size(dataTargetYog.X,1);
    dataTargetYog.y = dataMerged(:,end);
    dataTargetYog.max_x = dataSource.max_x;
    dataTargetYog.min_x = dataSource.min_x;
    dataTargetYog.mean_y = dataSource.mean_y;
    dataTargetYog.var_y = dataSource.var_y;
    
    gpTarget2_ytma = buildGP(dataTargetYog,paramTarget,boptions);

% recommend new sample 
    xnew2_ytma = recommendSampleTransfer(gpTarget2_ytma);
      
    newdata = newdataPoint(gpTarget2_ytma,testMode,xnew2_ytma);

    dataTargetYog.Xt = [dataTargetYog.Xt;newdata.X];
    dataTargetYog.y0 = [dataTargetYog.y0;newdata.y];
    
    [dataTargetYog.mean,dataTargetYog.var] = find_mean_var(dataTargetYog.y0,'target');
    dataTargetYog.y_re_computed = compute_response(dataTargetYog.y_dummy,dataTargetYog.mean,dataTargetYog.var,'target');
    response_target = compute_response(newdata.y,dataTargetYog.mean,dataTargetYog.var,'target');
    response_target(isnan(response_target)) = 0;
    dataTargetYog_new(:,end) = dataTargetYog.y_re_computed;
    dataTargetYog_new = [dataTargetYog_new;[newdata.X response_target]];
    
    dataTargetYog.y_dummy = [dataTargetYog.y_dummy;newdata.y];
   