function [ypmean, ypvar] = predictGP(gp,X1)

% msrSigma2 = gp.msrSigma2;
invK = gp.invK;
X = gp.X;
% X = gp.Xorg;
y = gp.y;
kmat = kernel(X,X1,gp.param);
ktinvK = kmat'*invK;
ypmean = ktinvK*y;
ypvar = (gp.param.kernelVar) - sum(ktinvK.*kmat',2);%(1) - sum(ktinvK.*kmat',2);%;

if ~isempty(find(ypvar < 0)) || sum(isnan(ypvar))
    ypvar(ypvar<0) = 0;
% %     keyboard
end