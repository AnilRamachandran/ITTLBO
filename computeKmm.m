function [ ret ] = computeKmm(Xbar, l, sigma, sigma0)

	m = size(Xbar, 1);
	Xbar = Xbar .* repmat(sqrt(l'), m, 1);
	Qbar = repmat(sum(Xbar.^2, 2), 1, m);
	distance = Qbar + Qbar' - 2 * Xbar * Xbar';
	ret = sigma * exp(-0.5 * distance) + eye(m) .* (sigma0 + sigma * 1e-10);
