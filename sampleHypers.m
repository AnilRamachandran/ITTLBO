
function [ l, sigma, sigma0 ] = sampleHypers(Xsamples, Ysamples, nSamples)

	% We specify a gamma prior for the noise

	pNoise = prior_logunif();

	% We specify the likelihood

	lik = lik_gaussian('sigma2_prior', pNoise);

	% We specify a gamma prior for the lengthscales

	meanG = 0.5;
	varianceG = 0.1;
	alpha = meanG^2 / varianceG;
	beta = meanG / varianceG;
	pLengthscale = prior_gamma('sh', alpha, 'is', beta);

	% We specify a gamma prior for the amplitud

	pAmplitud = prior_logunif();

	% We specify the covariance function

	gpcf = gpcf_sexp('lengthScale', 0.3 * ones(1, size(Xsamples, 2)), ...
		'lengthScale_prior', pLengthscale, 'magnSigma2_prior', pAmplitud);

	% We infer the hyper-parameters

	gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-10);
    
   
	% We sample from the hyper-parameters
    opt = optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
    gp = gp_optim(gp,Xsamples,Ysamples,'opt',opt);


    sigma0 = abs(normrnd(gp.lik.sigma2,0.001,nSamples,1));
    l_temp = (1 ./ gp.cf{ 1 }.lengthScale.^2);
    l = [];
    for ii = 1:size(l_temp,2)
        each_l = abs(normrnd(l_temp(ii),0.001,nSamples,1));
        l = [l each_l];
    end
    sigma = abs(normrnd(gp.cf{ 1 }.magnSigma2,0.001,nSamples,1));

