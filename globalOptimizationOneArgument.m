function [ optimum ] = globalOptimizationOneArgument(target, xmin, xmax, guesses)

	d = size(xmin, 1);

	% We evaluate the objective in a grid and pick the best location

	gridSize = 1000;

	Xgrid = repmat(xmin', gridSize, 1) + repmat((xmax - xmin)', gridSize, 1) .* rand(gridSize, d);
	Xgrid = [ Xgrid ; guesses ];

	y = target(Xgrid);
%     Problem.f = target;

	[ minValue minIndex ] = min(y);

	start = Xgrid(minIndex,:);

	% We optimize starting at the best location
    try
        optimum = fmincon(target, start, [], [], [], [], xmin, xmax, [], ...
                optimset('MaxFunEvals', 10, 'TolX', eps, 'Display', 'off', 'GradObj', 'on'));
    catch
        flag = 1;
    end

end
