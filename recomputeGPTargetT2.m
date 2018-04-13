function gpTarget = recomputeGPTargetT2(gpSource,gpTarget)

gpSource(1).alpha0 = 3;gpSource(1).beta0 = 3;
gpSource(2).alpha0 = 2;gpSource(2).beta0 = 1.5;
gpSource(3).alpha0 = 1;gpSource(3).beta0 = 0.2;
gpSource(4).alpha0 = 2.5;gpSource(4).beta0 = 2.5;
% gpSource(1).alpha0 = 2;gpSource(1).beta0 = 5;
% gpSource(2).alpha0 = 2;gpSource(2).beta0 = 5;
% gpSource(3).alpha0 = 2;gpSource(3).beta0 = 5;
% gpSource(4).alpha0 = 2;gpSource(4).beta0 = 5;

msrSigma2_source = [];
for ii = 1:length(gpSource)
    Ns = gpTarget.Ns;
    X = gpTarget.X(Ns+1:end,:);
    ys = predictGP(gpSource(ii),X);
    yt = gpTarget.y(Ns+1:end);
    Nt = length(yt);
    
    sum2 = (ys-yt)'*(ys-yt);
    alphan = gpSource(ii).alpha0 + (Nt/2);
    betan = gpSource(ii).beta0 + sum2/2;
    
    msrSigma2_sourcetemp(ii) = betan/(alphan+1);
    msrSigma2_source = [msrSigma2_source ones(1,size(gpSource(ii).X,1))*msrSigma2_sourcetemp(ii)];
end
sigma2Noise = diag([msrSigma2_source zeros(1,Nt)]);
gpTarget.invK = (gpTarget.K + sigma2Noise + diag(gpTarget.msrSigma2))\eye(size(gpTarget.K));%inv(gpTarget.K + sigma2Noise + diag(gpTarget.msrSigma2));