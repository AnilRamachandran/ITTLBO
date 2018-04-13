function [ ret ] = computeKnm(X, Xbar, l, sigma)

	n = size(X, 1);
	m = size(Xbar, 1);

	X = X .* repmat(sqrt(l'), n, 1);
	Xbar = Xbar .* repmat(sqrt(l'), m, 1);

	Q = repmat(sum(X.^2, 2), 1, m);
	Qbar = repmat(sum(Xbar.^2, 2)', n, 1);

	distance = Qbar + Q - 2 * X * Xbar';
	ret = sigma * exp(-0.5 * distance);
