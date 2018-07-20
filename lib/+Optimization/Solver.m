classdef (Abstract) Solver < handle

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties
		%	Define minimum and maximum number of iterations: the run loop will
		%	execute at least IterMin and no more than IterMax iterations
		IterMin(1,1) double {mustBePositive, mustBeNonNan} = 1;
		IterMax(1,1) double {mustBePositive, mustBeNonNan} = inf;

		%	The error function to be minimized: must accept a matrix of
		%	[PopLength x n] inputs and output [1 x n] error values
		ErrorFun(1,1) function_handle = @(p) sum(abs(p).^2);
		
		%	Function to constrain the population
		ConstrainFun(1,1) function_handle = @(p) p;

		%	Population range: the value of each gene must be greater than or
		%	equal to PopMin and less than or equal to PopMax. These values can
		%	be scalars or [PopLength x 1] arrays (none, one or both entities can
		%	be a scalar or an array)
		PopMin(:,1) double = 0;
		PopMax(:,1) double = 1;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (SetAccess=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess=private)
		%	The population array or gene values: matrix of size 
		%	[PopLength x PopNum]
		Population(:, :) double {mustBeFinite, mustBeNonempty} = 0;
		
		%	The number of individuals in the population
		PopNum(1,1) double {mustBeGreaterThanOrEqual(PopNum, 4), ...
			mustBeInteger} = 4;
		
		%	The length of (or number of genes in) an individual
		PopLength(1,1) double {mustBePositive, mustBeInteger} = 1
		
		%	The individual errors
		PopError(1, :) double;
		
		%	The optimal solution error
		SolnError(1,1) double;
		
		%	The optimal solution
		Solution(:, 1) double;
	
		%	Current iteration number
		Iter;
	end %	properties (SetAccess=private)
	
end