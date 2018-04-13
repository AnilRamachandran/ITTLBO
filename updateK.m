function gp = updateK(gp)

kvec = kernel(gp.X(1:end-1,:),gp.X(end,:),gp.param);
gp.K = [gp.K,kvec;kvec',gp.param.kernelVar];
gp.invK = pinv(gp.K);
% if cond(gp.K)>200
%     keyboard
% end