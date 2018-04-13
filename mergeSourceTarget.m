function [Xst,yst] = mergeSourceTarget(Xs,ys,Xt,ytn)
deltax = (0.25/2)*size(Xt,2);
yst = ys;
Xst = Xs;
M = size(Xs,1);
cnt = 1;
for ii = 1 : size(Xt,1)
    indexXt = find(sum(abs(Xs-kron(ones(size(Xs,1),1),Xt(ii,:))),2)<deltax);
    if isempty(indexXt)
        Xst(M + cnt,:) = Xt(ii,:);
        yst(M + cnt) = ytn(ii);
        cnt = cnt + 1;
    else
        yst(indexXt) = ytn(ii);
    end
end