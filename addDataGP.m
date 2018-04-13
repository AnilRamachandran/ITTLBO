function gp = addDataGP(gp, newdata)

xnew = newdata.X;
ynew = newdata.y;
gp.X = [gp.X;xnew];
gp.y = [gp.y;ynew];
gp.N = gp.N+1;
gp.msrSigma2 = [gp.msrSigma2;gp.msrSigma2scalar];
gp = updateK(gp);

if(isfield(gp.boptions,'SearchGrid'))
    SearchGrid = gp.boptions.SearchGrid;
    [SearchGrid] = removeXfromGrid(SearchGrid,newdata.X);
    gp.boptions.SearchGrid = SearchGrid;
end
