classdef Simplex < handle
	
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
		
		WeightedCentroid(1, 1) logical = false;
		Gamma(1, 1) double = 1;
		Vectorized(1, 1) logical = true;

		%	Scale factors
		Reflection(1,1) double {mustBePositive, mustBeFinite} = 1;
		Expansion(1,1) double {mustBeGreaterThan(Expansion, 1), ...
			mustBeFinite} = 2;
		Contraction(1, 1) double {mustBeLessThan(Contraction, 1), ...
			mustBeFinite} = .5;
		Shrink(1, 1) double {mustBeLessThan(Shrink, 1), ...
			mustBeFinite} = .5;
		
	
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

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (Hidden, Access=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess=private)
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		
		%	====================================================================
		%		C O N S T R U C T O R
		%	====================================================================
		function obj = Simplex(pop, fun, varargin)
			
			%	Initial conditions
			%	----------------------------------------------------------------
			if isempty(pop)
				%	Generate population using regular simplex
				obj.PopLength = FindArg(varargin, 'PopLength');
				pop = obj.RegularSimplex(obj.PopLength);
				
			elseif isvector(pop) && length(pop)>2
				obj.PopLength = length(pop);
				pop = pop + obj.RegularSimplex(obj.PopLength);
			else
				if size(pop, 2)~=(size(pop, 1)+1)
					error('Populations must be vector of size [N, N+1]');
				end
				
				%	Use defined initial population and extract population length
				%	and number
				obj.PopLength = size(pop, 1);
			end
			obj.PopNum = obj.PopLength + 1;
			
			%	Objective function
			obj.ErrorFun = fun;

			%	Get optional properties
			%	----------------------------------------------------------------
			obj.IterMin = FindArg(varargin, 'IterMin', obj.IterMin);
			obj.IterMax = FindArg(varargin, 'IterMax', obj.IterMax);

			obj.Vectorized = FindArg(varargin, 'Vectorized', obj.Vectorized);
			obj.ConstrainFun = FindArg(varargin, 'ConstrainFunction', ...
				obj.ConstrainFun);
			
			%	Check if weight centroid
			obj.WeightedCentroid = FindArg(varargin, ...
				'WeightedCentroid', obj.WeightedCentroid);
			obj.Gamma = FindArg(varargin, 'Gamma', obj.Gamma);
						
			%	Check for scaling parameter
			obj.Reflection = FindArg(varargin, 'Reflection', obj.Reflection);
			obj.Expansion = FindArg(varargin, 'Expansion', obj.Expansion);
			obj.Contraction = FindArg(varargin, 'Contraction', obj.Contraction);
			obj.Shrink = FindArg(varargin, 'Shrink', obj.Shrink);
				
			
			%	Check points are within bounds
			pop = obj.ConstrainFun(Coerce(pop, obj.PopMin, obj.PopMax));
					
			err = FindArg(varargin, 'PopError', []);
			if isempty(err)
				%	Calculate costs and sort
				err = obj.Evaluate(pop);
			end
			
			%	Sort population
			[err, ind] = sort(err);

			%	Update values
			obj.PopError = err;
			obj.Population = pop(:, ind);
			obj.Iter = 0;
			
			obj.SolnError = obj.PopError(1);
			obj.Solution = obj.Population(:, 1);
			
		end
				
		%	====================================================================
		%		M O V E   P O I N T
		%	====================================================================
		function Iterate(obj, Shift)
			
			%	Function to apply shift to cost function, default = no shift
			%	This can be used to implement simulated annealing
			if nargin<2 || isempty(Shift)
				Shift = @(F, s) F;
			end
			
			obj.Iter = obj.Iter + 1;
			
			%	Extract list of simplex vertices
			P0 = obj.Population;
			
			%	Extract cost values at each simplex vertex
			F0 = obj.PopError;
			
			%	Thermally shift each cost function
			G0 = Shift(F0, 1);
			
			%	Select current point
			P1 = P0(:, end);

			%	Calculate centroid of remaining points
			if obj.WeightedCentroid
				W = exp(-obj.Gamma*(F0-F0(1))./diff(F0([1 end])));
				Pcen = sum(P0(:, 1:end-1).*Row(W(1:end-1)), 2) ./ sum(W);
			else
				Pcen = mean(P0(:, 1:end-1), 2);
			end
			
			%	Reflect current point through centroid
			Vec = P1 - Pcen;
			Pr = obj.ConstrainFun(Coerce(Pcen - obj.Reflection.*Vec, ...
				obj.PopMin, obj.PopMax));
			
			%	Calculate cost
			Fr = obj.Evaluate(Pr);
			
			%	Subtract thermal shift
			Gr = Shift(Fr, -1);
			
			
			%	Check if current point better than second worst
			if Gr<G0(end-1)
				
				%	Check if current point is best
				if Gr<G0(1)
					
					%	Current point best --> Expand point
					Pe = obj.ConstrainFun(Coerce(Pcen ...
						+ obj.Expansion .* (Pr-Pcen), obj.PopMin, obj.PopMax));
					Fe = obj.Evaluate(Pe);
					Ge = Shift(Fe, -1);
					
					%	Check if expansion better than reflection
					if Ge<Gr
						%	Insert expanded point to beginning
						obj.Population = [Pe P0(:, 1:end-1)];
						obj.PopError = [Fe F0(1:end-1)];
					else
						%	Insert reflected point to beginning
						obj.Population = [Pr P0(:, 1:end-1)];
						obj.PopError = [Fr F0(1:end-1)];
					end %	if Ge<Gr
					
				else
					
					%	Insert reflected point into list
					ind = find(Gr<=G0, 1, 'first');
					if isempty(ind)
						ind = obj.PopNum;
					end
					obj.Population = [P0(:, 1:ind-1) Pr P0(:, ind:end-1)];
					obj.PopError = [F0(1:ind-1) Fr F0(ind:end-1)];
				end %	if
				
			else
				
				%	Reflected point not better: contract point
				Pc = obj.ConstrainFun(Coerce(Pcen ...
					+ obj.Contraction .* Vec, obj.PopMin, obj.PopMax));
				Fc = obj.Evaluate(Pc);
				Gc = Shift(Fc, -1);

				%	Check if contracted point is better than worst
				if Gc<G0(end)

					%	Insert contracted point into list
					ind = find(Gc<=G0, 1, 'first');
					if isempty(ind)
						ind = obj.PopNum;
					end
					obj.Population = [P0(:, 1:ind-1) Pc P0(:, ind:end-1)];
					obj.PopError = [F0(1:ind-1) Fc F0(ind:end-1)];

				else

					%	Shrink all points towards best
					Vec = P0(:, 1) - P0(:, 2:end);
					P1 = obj.ConstrainFun(Coerce(...
						P0(:, 1) + obj.Shrink * Vec, ...
						obj.PopMin, obj.PopMax));
					F1 = obj.Evaluate(P1);
					P0 = [P0(:, 1) P1];
					[obj.PopError, ind] = sort([F0(1) F1]);
					obj.Population = P0(:, ind);

				end %	Gc<G0(end)
			end %	if Gr<G0(end-1)
		
% 			%	Expand points along basis set
% 			sc = obj.Scale.Basis;
% 			if ~isnumeric(sc) || sc
% 				if ~isnumeric(sc)
% 					sc = sc();
% 				end
% 				if obj.WeightedCentroid
% 					W = exp(-obj.Gamma*(obj.PopError-obj.PopError(1)) ...
% 						./diff(obj.PopError([1 end])));
% 					Pcen = sum(obj.Population.*W, 2)./sum(W);
% 				else
% 					Pcen = mean(obj.Population, 2);
% 				end
% 				dP2 = sum(abs(obj.Population - Pcen).^2);
% 				RMS = sqrt(mean(dP2));
% 				obj.Population = obj.Coerce(obj.ConstrainFun(...
% 					obj.Population + sc*RMS * obj.Basis));
% 			end %	if ~isnumeric(sc) || sc

			if obj.PopError(1) < obj.SolnError
				obj.SolnError = obj.PopError(1);
				obj.Solution = obj.Population(:, 1);
			end
			
		end %	function
		
	
		%	====================================================================
		%		E V A L U A T E
		%	====================================================================
		function err = Evaluate(obj, points)
			if obj.Vectorized
				err = obj.ErrorFun(points);
			else
				err = zeros(1, obj.PopLength);
				for n=1:obj.PopLength
					err(n) = obj.ErrorFun(points(:, n));
				end
			end
		end
						
	end %	methods
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (Access=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access=private)
		
		%	====================================================================
		
		%	====================================================================

	end %	methods (Access=private)
	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (Static)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Static)
		
		%	====================================================================
		%		R E G U L A R   S I M P L E X   ( N )
		%	====================================================================
		%	Returns a regular simplex centred at (0, 0, ..., 0) and the first
		%	vertex located at (1, 0, ..., 0) with subsequent vertexes forming
		%	an angle acos -1/N with respect to each other.
		function v = RegularSimplex(N)
			
			%	Initialize vectors
			v = zeros(N, N+1);
			
			%	Set first vector
			v(1) = 1;
			
			%	Loop over subsequent vectors
			for n=2:N
			
				%	Apply dot product rule
				v(n-1, n:end) = -...
					(1/N + sum(v(1:n-2, n-1).*v(1:n-2, n), 1)) ...
					/ v(n-1, n-1);
				
				%	Apply Pythagoras with positive root
				v(n, n) = sqrt(1 - sum(v(1:n-1, n).^2, 1));
			end
			
			%	Last vertex uses negative root
			v(end, end) = -v(end, end-1);
	
		end
		
	end

end