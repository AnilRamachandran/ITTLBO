function [data] = newdataPoint(gp,testMode,xnew)

max_x = gp.max_x;
min_x = gp.min_x;
mean_y = gp.mean_y;
var_y = gp.var_y;

if testMode.flag    
    [xnew, y] = targetData(testMode,xnew);
    data.X = xnew;
    data.y0 = y;
    y = (y-mean_y)/sqrt(var_y);
    data.y = y;
end





