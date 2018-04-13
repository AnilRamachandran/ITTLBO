function [ optimum, minIndex] = globalOptimizationOneArgument(target, targetData_temp, guesses)

	

	% We evaluate the objective in a grid and pick the best location


	Xgrid = targetData_temp(:,1:3);
% 	Xgrid = [ Xgrid ; guesses ];

	y = target(Xgrid);
%     Problem.f = target;

	[ minValue minIndex ] = min(y);
    optimum = Xgrid(minIndex,:);
