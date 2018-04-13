function [fx, invK] = fit_rankGP(y, K, param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fx0 = rand(size(y));
if ~isrow(fx0)
    fx0 = fx0';
end

options = optimoptions('fminunc','GradObj','on','Display','off','MaxFunEvals',20,'MaxIter',20);%

dummy = 0;
[fx, fval] = fminunc(@(fx)mylikelihood(fx,y,K,param), fx0, options);

sq2ksigma = sqrt(2*param.kernelVar);

invK = zeros(size(K));

if ~isrow(y)
    y = y';
end
for ii = 1:length(y)
    z = (fx(ii)-fx)/sq2ksigma;
    invK(ii,:) = ((sign(y(ii)-y).*sign(y-y(ii))).*(((normpdf(z).^2)./(normcdf(z).^2)) + (z.*normpdf(z)./normcdf(z))))/(sq2ksigma^2);
end

invK = inv(K) + invK;

end

