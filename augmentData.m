function gp = augmentData(gp, newdata)

xnew = newdata.X;
ynew = newdata.y;
gp.X = [gp.X;xnew];
gp.y = [gp.y;ynew];
