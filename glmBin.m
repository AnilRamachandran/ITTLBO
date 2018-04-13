function [output] = glmBin(model,data,param,opt,mode)
%
% Generalized Linear Models for Binary outcomes
%
% Input:
%	data: 	each row is a data vector
%	label:  a column of 0/1 label
%	nIters: 	maximum number of iterations
%	l1_penalty: Laplacian penalty, often in the range [1e-5,1e-2];
%	l2_penalty: quadratic penalty, often in the range [1e-5,1e-2];
%	epsilon	:	threshold at which the learning stops if the relative improvement in log-likelihood falls below it.
%	report_interval: time interval to report the objective function (regularised likelihood). Set a very large number if you don't want to see report.
%
% Output:
%

[dataSize,N] = size(data);

if strcmp(upper(mode),'TRAIN') | strcmp(upper(mode),'TRAINING') | strcmp(upper(mode),'LEARN') | strcmp(upper(mode),'LEARNING')

	nIters			= opt.nIters;
	epsilon			= opt.epsilon;
	report_interval	= opt.report_interval;

	opt.model 	= model;
	opt.data 	= data;
	
	if isfield(opt,'isAddedBias') & opt.isAddedBias
		param 	= [0,param];
		opt.data = [ones(dataSize,1),data];
	end
	
	param = conjugate_gradient('glmBinGrad',param,opt,epsilon,nIters,report_interval);
	output = param;
else
	
	if length(param) == N + 1
		data = [ones(dataSize,1),data];
	end
	
	vals = data*param';

	if strcmp(upper(model),'PROBIT') | strcmp(upper(model),'GAUSS') | strcmp(upper(model),'GAUSSIAN') | strcmp(upper(model),'NORMAL')

		probs   = 1 - cdf('norm',-vals,0,1);

	elseif strcmp(upper(model),'GUMBEL')
	
		probs   = 1 - exp(-exp(vals));
		
	elseif strcmp(upper(model),'WEIBULL') | strcmp(upper(model),'GOMPERTZ')
	
		probs   = exp(-exp(vals));

	else %strcmp(upper(model),'LOGIT')

		probs	= 1./ (1 + exp(-vals));

	end

	output = probs;
end