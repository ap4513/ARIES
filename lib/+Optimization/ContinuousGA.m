classdef ContinuousGA < handle
	
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
% 		ErrorMin
% 		ErrorTol
% 

		%	Function to constrain the population
		ConstrainFun(1,1) function_handle = @(p) p;
		
		%	The number of parents to select for mating: the ParentNum
		%	individuals with the lowest error are selected based on rank
		%	weighting with replacement
		ParentNum(1, 1) double {mustBeGreaterThanOrEqual(ParentNum, 2), ...
			mustBeInteger} = 2;
		
		%	The number of children to generate each iteration: must be divisibe
		%	by 2 and less than PopNum-EliteNum. These children replace the
		%	ChildNum individuals with the highest error
		ChildNum(1, 1) double {mustBeGreaterThanOrEqual(ChildNum, 2), ...
			mustBeEven} = 2;
		
		%	The number of individuals passed to the next iteration without
		%	mutation or crossover: must be less than PopNum-ChildNum
		EliteNum(1,1) double {mustBeNonnegative, mustBeInteger} = 0;
		
		%	Mutation properties
		%		MutationRate is the probability of a particular gene mutating.
		%		If MutationAmplitude is infinte, then the mutated gene is given
		%		a random value between 0 & 1 (coerced if necessary). Otherwise
		%		the a normally distributed random offset with standard deviation
		%		MutationAmplitude is added to the gene. If MutationDecayRate is
		%		finite, the standard deviation reduces exponentially with
		%		iteration number.
		MutationRate(1,1) double {mustBeNonnegative} = inf;
		MutationAmplitude(:,1) double {mustBePositive, mustBeNonempty} = inf;
		MutationDecayRate(1,1) double {mustBePositive, mustBeNonNan} = inf;
		
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
		PopNum(1,1) double {mustBePositive, mustBeInteger} = 1;
		
		%	The length of (or number of genes in) an individual
		PopLength(1,1) double {mustBePositive, mustBeInteger} = 1
		
		%	The individual errors
		PopError(1, :) double;
		
		%	The optimal solution error
		SolnError(1,1) double;
		
		%	The optimal solution
		Solution(:, 1) double;
		
		%	Array of indices marking which individuals have been modified. Only
		%	these individuals have their error recalculated
		Modified(1, :) logical;
		
		%	Cumulative probability (rank weighting) for parent selection
		CumProb(1, :) double;
		
		%	Current iteration number
		Iter;
	end %	properties (SetAccess=private)
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (Dependent)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (Dependent)
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (SetAccess=private, Dependent)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess=private, Dependent)
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		
		%	====================================================================
		%		C O N T I N U O U S G A
		%	====================================================================
		%	Object constructor:
		%
		%		obj = ContinuousGA(pop, fun);
		%		obj = ContinuousGA([], fun, 'PopNum', n, 'PopLength', l);
		%		obj = ContinuousGA(___, PropName, PropValue, ...);
		%
		%	pop = Population array (PopLength x PopNum)
		%	fun = objective function to be minimised
		%	PropName, PropValue = List of name/value pairs
		%
		%	PropName			default	Descroption
		%	====================================================================
		%	IterMin				0		Minimum number of iterations
		%	IterMax				inf		Maximum number of iterations
		%	PopMin				0		Minimum population value
		%	PopMax				1		Maximum population value
		%	EliteNum			1		Number passed through without mutation
		%	ParentNum		   PopNum/2 Number of parents for mating
		%	ChildNum		   PopNum/2 Number of children produced by mating
		%	MutationRate		0.2		Probability of mutation
		%	MutationAmplitude	inf		Amplitude of mutation 
		%	MutationDecayRate	inf		Rate of decay of mutation amplitude
		%	ConstrainFunction	@(p) p	Constrains the population values
		function obj = ContinuousGA(pop, fun, varargin)
			
			%	Initial conditions
			%	----------------------------------------------------------------
			if isempty(pop)
				%	Generate random population from given population length and
				%	number
				obj.PopLength = FindArg(varargin, 'PopLength');
				obj.PopNum = FindArg(varargin, 'PopNum');
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
			obj.IterMin = FindArg(varargin, 'IterMin', obj.IterMin);
			obj.IterMax = FindArg(varargin, 'IterMax', obj.IterMax);
			
			obj.PopMin = FindArg(varargin, 'PopMin', obj.PopMin);
			obj.PopMax = FindArg(varargin, 'PopMax', obj.PopMax);
			obj.Population = Coerce(obj.Population, obj.PopMin, obj.PopMax);
			
% 			obj.ErrorMin = FindArg(varargin, 'ErrorMin', 0);
% 			obj.ErrorTol = FindArg(varargin, 'ErrorTol', 1e-6);
			
			obj.EliteNum = FindArg(varargin, 'EliteNum', 1);
			obj.ParentNum = FindArg(varargin, 'ParentNum', floor(obj.PopNum/2));
			obj.ChildNum = FindArg(varargin, 'ChildNum', 2*ceil(obj.PopNum/4));
			
			obj.MutationRate = FindArg(varargin, 'MutationRate', ...
				obj.MutationRate);
			obj.MutationDecayRate = FindArg(varargin, 'MutationDecayRate', ...
				obj.MutationDecayRate);
			obj.MutationAmplitude = FindArg(varargin, 'MutationAmplitude', ...
				obj.MutationAmplitude);
% 			obj.PlotOutput = FindArg(varargin, 'PlotOutput', false);
% 			obj.PlotRate = FindArg(varargin, 'PlotRate', 1);

			obj.ConstrainFun = FindArg(varargin, 'ConstrainFunction', ...
				obj.ConstrainFun);

			%	Reset GA
			obj.Reset;
			
		end %	function ContinuousGA
		
		%	====================================================================
		%		R E S E T
		%	====================================================================
		function Reset(obj)
			
			%	Set iteration counter to zero
			obj.Iter = 0;
			
			%	Reset optimal solution
			obj.SolnError = inf;
			
			%	Calculate errors from current population
			obj.Modified = true(1, obj.PopNum);
			obj.CalculateError;
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
			
			%	Perform mating ritual
			obj.Mate;
			
			%	Mutate population
			obj.Mutate;
			
			%	Constrain values
			obj.Population = obj.ConstrainFun(obj.Population);
			
			%	Update errors
			obj.CalculateError;
		end
		
	end %	methods
	
	%	========================================================================
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S   ( S e t / G e t )
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		
		%	====================================================================
		%		P A R E N T   N U M
		%	====================================================================
		function set.ParentNum(obj, val)
			val = Coerce(val, 2, obj.PopNum); %#ok<MCSUP>
			obj.ParentNum = val;
			obj.CumProb = cumsum((val-(0:val-1)) * 2/(val*(val+1))); %#ok<MCSUP>
		end
		
		%	====================================================================
		%		C H I L D   N U M
		%	====================================================================
		function set.ChildNum(obj, val)
			if val>obj.PopNum-obj.EliteNum %#ok<MCSUP>
				error('Number of children plus elite must not exceed population number');
			end
			obj.ChildNum = val;
		end
		
		%	====================================================================
		%		E L I T E   N U M
		%	====================================================================
		function set.EliteNum(obj, val)
			if val>obj.PopNum-obj.ChildNum %#ok<MCSUP>
				error('Number of children plus elite must not exceed population number');
			end
			obj.EliteNum = val;
		end
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S   (Access=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access=private)
				
		%	====================================================================
		%		C A L C U L A T E   E R R O R
		%	====================================================================
		function CalculateError(obj)
			%	Only calculate error for modified populations
			obj.PopError(obj.Modified) = ...
				obj.ErrorFun(obj.Population(:, obj.Modified));
			
			%	Sort population
			[obj.PopError, ind] = sort(obj.PopError, 'ascend');
			obj.Population = obj.Population(:, ind);
			
			%	Check for better solution
			if obj.PopError(1)<=obj.SolnError
				obj.SolnError = obj.PopError(1);
				obj.Solution = obj.Population(:, 1);
			end
			
			%	Reset modified flag (all population errors have been calculated)
			obj.Modified = false(1, obj.PopNum);
		end
		
		%	====================================================================
		%		M A T E
		%	====================================================================
		function Mate(obj)
			%	Initialize children array
% 			Children = zeros(obj.PopLength, obj.ChildNum);
			
			%	Loop through mating ritual
% 			for n=1:2:obj.ChildNum
% 				Children(:, n+(0:1)) = obj.CrossOver;
% 			end
			
			%	Update populations and modified indices
			obj.Population(:, end-obj.ChildNum+1:end) = Coerce(...
				obj.CrossOver, obj.PopMin, obj.PopMax);
			obj.Modified(end-obj.ChildNum+1:end) = true;
		end
		
		%	====================================================================
		%		G E T   P A R E N T   I N D E X
		%	====================================================================
		function ind = GetParentIndex(obj, n)
			ind = 1 + sum(rand(n, 1)>obj.CumProb, 2);
		end
		
		%	====================================================================
		%		C R O S S   O V E R
		%	====================================================================
		function Children = CrossOver(obj)
			%	Get indices of parents
			ind = reshape(obj.GetParentIndex(obj.ChildNum), [], 2);
			N = size(ind, 1);
			
% 			%	Ensure parents are unique
% 			while ind(1)==ind(2)
% 				ind(2) = obj.GetParentIndex(1);
% 			end
			
			%	Extract parents
			P1 = obj.Population(:, ind(:, 1));
			P2 = obj.Population(:, ind(:, 2));
			
			%	Index of crossover
% 			n = (1:obj.PopLength).';
			[~, n] = sort(rand(obj.PopLength, N));
% 			n0 = randi(obj.PopLength, [1 N]);
			n0 = obj.PopLength/2;
% 			nw = randi(obj.PopLength, [1 N]);
% 			nw = abs(randn(obj.PopLength, N)) * N/4;
			nw = abs(randn(1, N)) * N/4;
			
			beta = .5 + tanh(-((n-n0)./nw))/2;
			beta1 = 1 - beta;
			
			Children = [(P1.*beta + P2.*beta1) (P1.*beta1 + P2.*beta)];
			
% 			%	Random ratio at crossover point
% 			beta = rand;
% 			
% 			%	Generate Children
% 			Children = [ ...
% 				[P1(1:n-1); ...
% 				(1-beta)*P1(n) + beta*P2(n); ...
% 				P2(n+1:end)] ...
% 				...
% 				[P2(1:n-1); ...
% 				(1-beta)*P2(n) + beta*P1(n); ...
% 				P1(n+1:end)]];			
		end
		
		%	====================================================================
		%		M U T A T E
		%	====================================================================
		function Mutate(obj)
			
			%	Number of mutatable populations
			N = obj.PopNum-obj.EliteNum;
			
			%	Apply shift by an amount proportional to MR and in random
			%	direction
			if ~isfinite(obj.MutationRate)
				%	Random vector
				val = randn(obj.PopLength, N);
				
				%	Scale vector
				val = obj.MutationAmplitude.*(val./sqrt(sum(val.^2)));
				if isfinite(obj.MutationDecayRate)
					val = val .* exp(-obj.Iter/obj.MutationDecayRate);
				end
				
				%	Apply shift to population
				obj.Population(:, obj.EliteNum+1:end) = Coerce( ...
					obj.Population(:, obj.EliteNum+1:end) + val, ...
					obj.PopMin, obj.PopMax);
				
				%	All genes mutated
				obj.Modified(obj.EliteNum+1:end) = true;
				
				%	End of mutation
				return
			end
			
			%	Alternatively mutate individual genes
			%	Get indices of mutations
			ind = Column([false(obj.PopLength, obj.EliteNum) ...
				(rand(obj.PopLength, N)<obj.MutationRate)]);
			
			%	Total number of mutations
			M = sum(ind);
			
			%	Only bother if there is a mutation
			if M
				
				%	Need to get bounds of mutated object
				[ny, nx] = ind2sub([obj.PopLength, obj.PopNum], find(ind));
				if isscalar(obj.PopMin)
					pmin = obj.PopMin;
				else
					pmin = obj.PopMin(ny);
				end
				
				if isscalar(obj.PopMax)
					pmax = obj.PopMax;
				else
					pmax = obj.PopMax(ny);
				end
				
				%	If finite, then add a normally distributed random offset,
				%	otherwise mutate to uniformly distributed value
				if ~isscalar(obj.MutationAmplitude) ...
						|| isfinite(obj.MutationAmplitude)
					
					%	Random offset
					if isscalar(obj.MutationAmplitude)
						val = obj.MutationAmplitude;
					else
						val = obj.MutationAmplitude(ny);
					end
% 					val = val .* randn(M, 1) .* (1-exp(-(3*obj.SolnError).^.25));
					
					val = val .* randn(M, 1);
					
					%	Reduce offset with iteration if required
					if isfinite(obj.MutationDecayRate)
						val = val * exp(-(obj.Iter/obj.MutationDecayRate));
					end
					
					%	Add random offset & coerce to within range
					obj.Population(ind) = Coerce(obj.Population(ind) + val, ...
						pmin, pmax);
				else
					val = rand(M, 1);
					if isfinite(obj.MutationDecayRate)
						val = obj.Population(ind) + (2*val - 1) ...
							* exp(-(obj.Iter/obj.MutationDecayRate).^2);
					end
					%	Make a random value
					obj.Population(ind) = Coerce(val, pmin, pmax);
				end
				
				obj.Modified(nx) = true;
			end %	if M
		end %	function Mutate
				
		%	====================================================================
		%		P L O T
		%	====================================================================
		function Plot(obj) %#ok<MANU>
% 			ax = axis;
% 			x = linspace(-1, 1, obj.PopLength)';
% 			h = get(gca, 'Children');
% 			plot(x, (obj.Population), ':', ...
% 				h(1).XData, [spline(x, obj.Solution, h(1).XData(:)) h(1).YData(:)], '+-');
% % 			axis tight
% % 			axis([-1 1 -2 2])
% 			axis(ax);
% 			title([obj.Iter obj.PopNum toc obj.SolnError max(obj.RelError)]);
% 			drawnow;
		end %	function Plot
				
	end %	methods (Access=private)
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S   (Static)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Static)
	end
end %	classdef
