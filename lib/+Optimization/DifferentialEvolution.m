classdef DifferentialEvolution < handle

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

		%	The number of vectors to use for the mutation
		VecNum(1, 1) double {mustBePositive, mustBeFinite, mustBeInteger} = 1;
		
		%	Crossover rate
		CrossoverRate(1,1) double {mustBeFinite, ...
			mustBeGreaterThanOrEqual(CrossoverRate, 0), ...
			mustBeLessThanOrEqual(CrossoverRate, 1)} = .5;

		%	Mutation factor
		MutationFactor(1,1) double {mustBeFinite, mustBeNonnegative} = .5;
		
		%	Decay rate
		DecayRate(1, 1) double {mustBePositive} = inf;
		
		%	Population range: the value of each gene must be greater than or
		%	equal to PopMin and less than or equal to PopMax. These values can
		%	be scalars or [PopLength x 1] arrays (none, one or both entities can
		%	be a scalar or an array)
		PopMin(:,1) double = 0;
		PopMax(:,1) double = 1;
	end %	properties
	
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
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		
		%	====================================================================
		%		C O N S T R U C T O R
		%	====================================================================
		function obj = DifferentialEvolution(pop, fun, varargin)
			
			%	Initial conditions
			%	----------------------------------------------------------------
			if isempty(pop)
				%	Generate random population from given population length and
				%	number
				obj.PopLength = FindArg(varargin, 'PopLength');
				obj.PopNum = FindArg(varargin, 'PopNum', obj.PopNum);
				obj.Population = rand(obj.PopLength, obj.PopNum);
			else
				%	Use defined initial population and extract population length
				%	and number
				obj.Population = pop;
				[obj.PopLength, obj.PopNum] = size(obj.Population);
			end
			
			%	Objective function
			obj.ErrorFun = fun;
			
			%	Get optional properties
			%	----------------------------------------------------------------
			obj.ConstrainFun = FindArg(varargin, 'ConstrainFunction', ...
				obj.ConstrainFun);
			obj.IterMin = FindArg(varargin, 'IterMin', obj.IterMin);
			obj.IterMax = FindArg(varargin, 'IterMax', obj.IterMax);
			
			obj.PopMin = FindArg(varargin, 'PopMin', obj.PopMin);
			obj.PopMax = FindArg(varargin, 'PopMax', obj.PopMax);
			obj.Population = obj.ConstrainFun(...
				Coerce(obj.Population, obj.PopMin, obj.PopMax));
			
			obj.VecNum = FindArg(varargin, 'VecNum', obj.VecNum);
			
			obj.CrossoverRate = FindArg(varargin, 'CrossoverRate', ...
				obj.CrossoverRate);
			obj.MutationFactor = FindArg(varargin, 'MutationFactor', ...
				obj.MutationFactor);
			obj.DecayRate = FindArg(varargin, 'DecayRate', obj.DecayRate);


			%	Reset GA
			obj.Reset(varargin);
			
		end
		
		
		%	====================================================================
		%		R E S E T
		%	====================================================================
		function Reset(obj, varargin)
			
			%	Set iteration counter to zero
			obj.Iter = 0;
			
			%	Reset optimal solution
			obj.SolnError = inf;
			
			%	Calculate errors from current population
			err = FindArg(varargin, 'PopError', []);
			if isempty(err)
				%	Calculate costs and sort
				obj.PopError = obj.ErrorFun(obj.Population);
			else
				obj.PopError = err;
			end
			obj.Sort;

% 			obj.CalculateError;
		end %	function Reset

	
		%	====================================================================
		%		R U N
		%	====================================================================
		function Run(obj) %#ok<MANU>
			error('Not implemented');
			Converged = false; %#ok<UNRCH>
			
			while ~Converged
			end
						
		end %	function Run

	
		%	====================================================================
		%		I T E R A T E
		%	====================================================================
		%	Single iteration step
		function Iterate(obj)
			%	Increment iteration counter
			obj.Iter = obj.Iter + 1;
			
			if obj.MutationFactor>0
				
				%	Generate mutated vectors
				ind = obj.GetPopIndices;
% 				V = Coerce(obj.Population(:, ind(1, :)) + obj.MutationFactor ...
% 					* sum(obj.Population(:, ind(2:2:end, :)) ...
% 					- obj.Population(:, ind(3:2:end, :)), 2), ...
% 					obj.PopMin, obj.PopMax);
				V = Coerce(obj.Population(:, ind(1, :)) + obj.MutationFactor ...
					* squeeze(sum(diff(reshape( ...
					obj.Population(:, ind(2:end, :)), obj.PopLength, 2, ...
					obj.VecNum, obj.PopNum), [], 2), 3)), ...
					obj.PopMin, obj.PopMax);
				
				%	Perform crossover
				if obj.CrossoverRate>0
					ind = (rand(obj.PopLength, obj.PopNum)<obj.CrossoverRate);
					V(ind) = obj.Population(ind);
				end
				
				m = mean(V, 2);
				if isfinite(obj.DecayRate)
					V = Coerce(m + (V-m) ./ (1-exp(-(obj.Iter/obj.DecayRate)...
						.*(max(V, [], 2)-min(V, [], 2)).^2)), ...
						obj.PopMin, obj.PopMax);
				end
				
				%	Constrain values
				obj.Population = obj.ConstrainFun(obj.Population);

				%	Calculate errors
				err = obj.ErrorFun(V);

				%	Update populations
				ind = (err<obj.PopError);
				obj.Population(:, ind) = V(:, ind);
				obj.PopError(ind) = err(ind);
				obj.Sort;
			end
			
			%	Update errors
% 			obj.CalculateError;
		end
		
	end %	methods
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S   ( S e t / G e t )
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		
		%	====================================================================
		%		V E C   N U M
		%	====================================================================
		function set.VecNum(obj, val)
			obj.VecNum = Coerce(val, 1, floor(obj.PopNum/2 - 1)); %#ok<MCSUP>
		end
		
	end %	methods
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S   (Access=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access=private)
		
		%	====================================================================
		%		S O R T
		%	====================================================================
		%	Sort populations according to error and check if solution has been
		%	improved
		function Sort(obj)
			[obj.PopError, ind] = sort(obj.PopError, 'ascend');
			obj.Population = obj.Population(:, ind);
			
			if obj.PopError(1)<obj.SolnError
				obj.SolnError = obj.PopError(1);
				obj.Solution = obj.Population(:, 1);
			end
			
		end
		
		%	====================================================================
		%		G E T   P O P   I N D I C E S
		%	====================================================================
		%	Generate list of population indices for generating mutation vectors
		function ind = GetPopIndices(obj)
			%	Number of vector samples
			NV = obj.VecNum*2 + 1;
			
			%	Number of individuals
			NP = obj.PopNum;
			
			%	If number of individuals small, then use sort routine, otherwise
			%	use iteration routine - this ratio was found emperically and is
			%	very approximate, but good enough
			if NP<10*NV
				%	Get list of random indices
				[~, ind] = sort(rand(NP-1, NP));
				
				%	Select first set of indices
				ind = ind(1:NV, :);
				
				%	Ensure not equal to index number
				n = (ind>=1:NP);
				ind(n) = ind(n) + 1;
				
			else
				
				%	Generate random indices
				ind = randi(NP-1, [NV NP]);
				
				%	Index vector
				N = 1:NP;
				
				%	Loop through each vector index
				for nv=2:NV
					
					%	Extract current indices
					n = ind(nv, :);
					
					%	Find which indices
					m = (n==N) | any(n==ind(1:nv-1, :), 1);
					while any(m)
						n(m) = randi(NP, [1 sum(m)]);
						m = (n==N) | any(n==ind(1:nv-1, :), 1);
					end
					ind(nv, :) = n;			
				end
			end
			
		end
				
		%	====================================================================
		%		C A L C U L A T E   E R R O R
		%	====================================================================
% 		function CalculateError(obj, pop)
% 			%	Only calculate error for modified populations
% 			obj.PopError(obj.Modified) = ...
% 				obj.ErrorFun(obj.Population(:, obj.Modified));
% 			
% 			%	Sort population
% 			[obj.PopError, ind] = sort(obj.PopError, 'ascend');
% 			obj.Population = obj.Population(:, ind);
% 			
% 			%	Check for better solution
% 			if obj.PopError(1)<=obj.SolnError
% 				obj.SolnError = obj.PopError(1);
% 				obj.Solution = obj.Population(:, 1);
% 			end
% 			
% 			%	Reset modified flag (all population errors have been calculated)
% 			obj.Modified = false(1, obj.PopNum);
% 		end
		
	end %	methods (Access=private)
	
	
end %	classdef