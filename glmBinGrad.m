function [ll,dl] = glmBinGrad(param,opt)

model		= opt.model;
data		= opt.data;
label		= opt.label;
l1_penalty	= opt.l1_penalty;
l2_penalty	= opt.l2_penalty;


dataSize = length(label);

vals = data*param';

if strcmp(upper(model),'PROBIT') | strcmp(upper(model),'GAUSS') | strcmp(upper(model),'GAUSSIAN') | strcmp(upper(model),'NORMAL')

	probs   = 1 - cdf('norm',-vals,0,1);
	dl = ((label - probs).*  normpdf(-vals) ./ (1e-10 + probs .* (1-probs)))'*data;
	
elseif strcmp(upper(model),'GUMBEL')

	probs   = 1 - exp(-exp(vals));
	dl = ((label - probs).* exp(vals) ./ (1e-10 + probs))'*data;

elseif strcmp(upper(model),'WEIBULL') | strcmp(upper(model),'GOMPERTZ')

	probs   = exp(-exp(vals));
	dl = -((label - probs).* exp(vals) ./ (1e-10 + 1-probs))'*data;

else %LOGIT as default

	probs	= 1./ (1 + exp(-vals));
	dl = (label - probs)'*data;

end

ll = label'*log(1e-10 + probs) + (1-label)'*log(1e-10 + 1-probs);

epsilon = 1e-10;

ll = ll/dataSize - 0.5*l2_penalty*param*param' - l1_penalty*sum(sqrt(param.^2 + epsilon));
dl = dl ./ dataSize - l2_penalty*param - l1_penalty*param ./ sqrt(param.^2 + epsilon);
