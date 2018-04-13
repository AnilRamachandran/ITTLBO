function [llk, llkg] = mylikelihood(fx,y,K,param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fx = fx';
invK = inv(K);
llk = 0.5*fx'*invK*fx;
llkg = zeros(1,length(y));
val1 = 0.0;
sq2ksigma = sqrt(2*param.kernelVar);
Ny = length(y);

for ii = 1:Ny
    z = (fx(ii)-fx(ii+1:Ny))/sq2ksigma;
    val1 = val1 + sum(log(normcdf(z)));
end


llk = -val1 + llk;

val2 = 0.0;
for ii = 1:Ny
 
    z = (fx(ii)-fx(ii+1:Ny))/sq2ksigma;
    sk = -1*sign(y(ii)-y(ii+1:Ny));
    if (~isrow(sk))
        sk = sk';
    end
    val2 = val2 + (sk*(normpdf(z)./normcdf(z)))/sq2ksigma;

    llkg(ii) = val2;
end 


llkg = llkg + (invK*fx)';
dummy = 0;
