classdef SimulatedAnnealing < uiobjects.Object & handle
	
	%	========================================================================
	%	P R O P E R T I E S
	%	========================================================================
	
	properties (Dependent, SetAccess=private)
		Solution double;
		OptimalSolution double;
		Temperature double;
		MinTemperature double;
		Iteration double;
		Cost double;
		OptimalCost double;
		StepSize double;
		StepSizeScaling double;
		Index double;
		UseStepSize logical;
		IterMax;
		IterMin;
	end
	
	properties (Hidden)
		Soln Optimization.Parameter;
		StpSz Optimization.Parameter;
		PlotRate double						%	Rate at which data is plotted
		PlotData logical;					%	Perform plotting?
		RunParallel logical;

 		%	Functions
		CostFn function_handle;				%	Function to minimize
		AnnealFn function_handle;			%	Temperature cooling function
		IndexFn function_handle;			%	Which param to change
		NeighbourFn function_handle;		%	Next trial solution generator
		AcceptanceFn function_handle;		%	Probability of acceptance
		ConvergenceFn function_handle;		%	Determine if converged
		PlotFn function_handle;				%	Plot data

	end %	properties
	
	%	========================================================================
	
	properties (Hidden)
% 		SolutionAcceptance double;			%	Number of solutions accepted
% 		TempIter double;					%	Iterations per temperature cycle
	end
	
	%	========================================================================

	properties (SetAccess=private, Hidden)
		OptSoln Optimization.Parameter;
		Cst Optimization.Parameter;
		Temp Optimization.Parameter;
		Iter Optimization.Parameter;
		Indx Optimization.Parameter;		%	Current solution vector index to change				
	end %	properties (SetAccess=private)
	
	
	
	%	========================================================================
	%	M E T H O D S
	%	========================================================================

	methods
		
		%	====================================================================

		function val = get.Solution(obj)
			val = obj.Soln.Value;
		end
		
		%	--------------------------------------------------------------------
		
		function val = get.OptimalSolution(obj)
			val = obj.OptSoln.Value;
		end
		
		%	--------------------------------------------------------------------
		
		function val = get.Temperature(obj)
			val = obj.Temp.Value;
		end
		
		%	--------------------------------------------------------------------
		
		function val = get.Iteration(obj)
			val = obj.Iter.Value;
		end
		
		%	--------------------------------------------------------------------

		function val = get.Cost(obj)
			val = obj.Cst.Value;
		end
		
		%	--------------------------------------------------------------------

		function val = get.OptimalCost(obj)
			val = obj.Cst.Lowest;
		end
		
		%	--------------------------------------------------------------------
		
		function val = get.StepSize(obj)
			val = obj.StpSz.Value;
		end

		%	--------------------------------------------------------------------
		
		function val = get.StepSizeScaling(obj)
			val = obj.StpSz.GetUserData('Scaling');
		end

		%	--------------------------------------------------------------------
		
		function val = get.Index(obj)
			val = obj.Indx.Value;
		end
		
		%	--------------------------------------------------------------------

		function val = get.UseStepSize(obj)
			val = obj.StpSz.GetUserData('Utilize');
		end
		
		%	--------------------------------------------------------------------
		
		function val = get.IterMax(obj)
			val = obj.Iter.Maximum;
		end
		
		%	--------------------------------------------------------------------
		
		function val = get.IterMin(obj)
			val = obj.Iter.Minimum;
		end
		
		%	--------------------------------------------------------------------

		function val = get.MinTemperature(obj)
			val = obj.Temp.Minimum;
		end
		
		%	====================================================================

		function set.Solution(obj, val)
			obj.Soln.SetValue(val);
		end
		
		%	--------------------------------------------------------------------
		
		function set.OptimalSolution(obj, val)
			obj.OptSoln.SetValue(val);
		end
		
		%	--------------------------------------------------------------------
		
		function set.Temperature(obj, val)
			obj.Temp.SetValue(val);
		end
		
		%	--------------------------------------------------------------------
		
		function set.Iteration(obj, val)
			obj.Iter.SetValue(val);
		end
		
		%	--------------------------------------------------------------------

		function set.Cost(obj, val)
			obj.Cst.SetValue(val);
		end
		
		%	--------------------------------------------------------------------
		
		function set.StepSize(obj, val)
			obj.StpSz.SetValue(val);
		end

		%	--------------------------------------------------------------------

		function set.StepSizeScaling(obj, val)
			obj.StpSz.SetUserData('Scaling', val);
		end
		
		%	--------------------------------------------------------------------
		
		function set.Index(obj, val)
			obj.Indx.SetValue(val);
		end
		%	====================================================================

		function obj = SimulatedAnnealing(soln, costfn, varargin)
			
			
			%	Set current solution
			%	----------------------------------------------------------------
			Opt = obj.SetOptions( struct( ...
					'Maximum', inf, ...
					'Minimum', -inf, ...
					'SaveHistory', false, ...
					'SaveHighest', false, ...
					'SaveLowest', false, ...
					'NotifyChange', false, ...
					'CoerceValues', true, ...
					'UserData', {[]}), ...
				obj.FindArg(varargin, 'SolutionOptions'));
			obj.Soln = Optimization.Parameter(soln, Opt{:});
			
			
			%	Set optimal solution (equal to current)
			%	----------------------------------------------------------------
			Opt = obj.SetOptions( struct( ...
					'SaveHistory', false, ...
					'SaveHighest', false, ...
					'SaveLowest', false, ...
					'NotifyChange', true), ...
				obj.FindArg(varargin, 'OptimalSolutionOptions'));
			obj.OptSoln = Optimization.Parameter(obj.Soln, ...
				'CoerceValues', false, Opt{:});
			
			
			%	Set cost function and current cost
			%	----------------------------------------------------------------
			obj.CostFn = costfn;
			Opt = obj.SetOptions( struct( ...
					'SaveHistory', false, ...
					'NotifyChange', false, ...
					'UserData', {[]}), ...
				obj.FindArg(varargin, 'CostOptions'));
			obj.Cst = Optimization.Parameter(costfn(soln), ...
				'CoerceValues', false, 'SaveLowest', true, Opt{:});
			
			
			%	Set up iteration counter
			%	----------------------------------------------------------------
			Opt = obj.SetOptions( struct( ...
					'Maximum', inf, ...
					'Minimum', 1, ...
					'SaveHistory', false, ...
					'NotifyChange', false), ...
				obj.FindArg(varargin, 'IterationOptions'));
			obj.Iter = Optimization.Parameter(0, ...
				'CoerceValues', false, Opt{:});
			
			
			%	Set up temperature
			%	----------------------------------------------------------------
			Opt = obj.SetOptions( struct( ...
					'Maximum', inf, ...
					'Minimum', 0, ...
					'SaveHistory', true, ...
					'NotifyChange', true, ...
					'Initial', 1), ...
				obj.FindArg(varargin, 'TemperatureOptions'));
			obj.Temp = Optimization.Parameter('CoerceValues', false, Opt{:});
	
			
			%	Set step size
			%	----------------------------------------------------------------
			Opt = obj.FindArg(varargin, 'StepSizeOptions');
			USS = obj.FindArg(Opt, 'Utilize', true);
			SSS = obj.FindArg(Opt, 'Scaling', .3);
			if SSS>1-1e-3
				SSS = 1-1e-3;
			elseif SSS<1e-3
				SSS = 1e-3;
			end
			Opt = obj.SetOptions( struct( ...
					'Maximum', inf, ...
					'Minimum', -inf, ...
					'SaveHistory', false, ...
					'SaveHighest', false, ...
					'SaveLowest', false, ...
					'NotifyChange', false, ...
					'CoerceValues', false, ...
					'Initial', ones(length(obj.Soln), 1), ...
					'UserData', {{'Utilize', USS, 'Scaling', SSS}}), ...
				Opt);
			obj.StpSz = Optimization.Parameter(Opt{:});
			if length(obj.StpSz)~=length(obj.Soln)
				obj.StepSize = repmat(obj.StepSize, size(obj.Soln));
			end
			
			%	Set index
			%	----------------------------------------------------------------
			Opt = obj.SetOptions( struct( ...
					'Maximum', length(obj.Soln), ...
					'Minimum', 1, ...
					'SaveHistory', false, ...
					'Initial', 1), ...
				obj.FindArg(varargin, 'IndexOptions'));
			obj.Indx = Optimization.Parameter(Opt{:});
		
			%	Extract additional parameters
			%	----------------------------------------------------------------
			obj.PlotRate = FindArg(varargin, 'PlotRate', 1);
			obj.PlotData = FindArg(varargin, 'PlotData', false);
			obj.RunParallel = FindArg(varargin, 'RunParallel', false);
			
			
			%	Set function handles
			%	----------------------------------------------------------------
			obj.AnnealFn = FindArg(varargin, ...
				'AnnealFunction', ...
				@(obj) obj.DefaultAnnealingScheduleFn);
			
			obj.NeighbourFn = FindArg(varargin, ...
				'NeighbourFunction', @(obj) obj.DefaultNeighbourFn);
			
			obj.IndexFn = FindArg(varargin, ...
				'IndexFunction', @(obj) obj.DefaultIndexFn);
			
			obj.AcceptanceFn = FindArg(varargin, ...
				'AcceptanceFunction', ...
				@(obj, cost)obj.DefaultAcceptanceFn(cost));

			obj.ConvergenceFn = FindArg(varargin, ...
				'ConvergenceFunction', @(obj) obj.DefaultConvergenceFn);
			
			obj.PlotFn = FindArg(varargin, ...
				'PlotFunction', @(obj) obj.DefaultPlotFn);

% 			obj.Reset;

		end

		%	====================================================================
		
		function Run(obj)
			stop = false;
			while ~stop
			
				obj.Iteration = obj.Iteration + 1;

				%	Get index
				ind = obj.IndexFn(obj);
				obj.Index = ind;

				%	Update solution
				sol = obj.NeighbourFn(obj);

				%	Calculate cost
				cost = obj.CostFn(sol);

				if obj.AcceptanceFn(obj, cost)

					%	Keep new solution
					%	----------------------------------------------------
					obj.Solution = sol;
					obj.Cost = cost;

					%	Update optimal solution
					if cost<=obj.Cst.Lowest
						obj.OptimalSolution = sol;
					end

					%	Update step size
					%		Increase step size in current direction
					obj.StepSize(ind) = obj.StepSize(ind)*(1+obj.StepSizeScaling);
					
				else
					%	Keep old solution
					%	----------------------------------------------------
					
					%	Update step size
					%		Decrease step size and randomize direction
					obj.StepSize(ind) = obj.StepSize(ind) ...
						* (-1)^randi(2)*obj.StepSizeScaling + eps;
					obj.Cost = obj.Cost;

				end

				if obj.PlotData && ~mod(obj.Iteration, obj.PlotRate)
					obj.PlotFn(obj);
				end

				stop = obj.ConvergenceFn(obj);
				
				obj.Temperature = obj.AnnealFn(obj);
				
			end
		end
		
		%	====================================================================

		function Reset(obj)
% 			obj.Temperature = obj.StartTemperature;
% 			obj.OptimalSolution = obj.Solution;
% 			obj.OptimalCost = obj.Cost;
% 			obj.SolutionCost = repmat(obj.Cost, [obj.SolutionLength 1]);
% 			obj.PreviousSolutionCost = zeros(obj.SolutionLength, 1);
% 			obj.Iter = 0;
% 			obj.Index = 1;
% 			obj.SolutionStep = ones(obj.SolutionLength, 1);
% 			obj.SolutionAcceptance = 0.5*obj.StepCycles*obj.SolutionStep;
% 			
% 			if obj.SaveCostHistory
% 				obj.CostHistory = repmat(obj.Cost, [100 1]);
% 			end
		end
		
		%	====================================================================
		
		function SetOptimalSolution(obj, val)
% 			val = Column(val);
% 			if length(val)~=obj.SolutionLength
% 				error('Does not match current size');
% 			end
% 			
% 			cost = obj.CostFn(val);
% 			if cost<=obj.OptimalCost
% 				obj.OptimalCost = cost;
% 				obj.OptimalSolution = val;
% 			end
		end

		%	====================================================================

		function SetCurrentSolution(obj, val)
% 			val = Column(val);
% 			if length(val)~=obj.SolutionLength
% 				error('Does not match current size');
% 			end
% 			
% 			cost = obj.CostFn(val);
% 			obj.Cost = cost;
% 			obj.Solution = val;
		end
		
		%	====================================================================
		
	end %	methods
	
	
	
	
	%	========================================================================
	%	M E T H O D S (Access=public)
	%	========================================================================
		
	methods (Access=public)
		
		%	====================================================================

		function Options = SetOptions(obj, default, varargin)
			fn = fieldnames(default);
			val = struct2cell(default);
			Options = cell(1, 2*length(fn));
			for n=1:length(fn)
				Options{2*n-1} = fn{n};
				Options{2*n} = obj.FindArg(varargin{:}, fn{n}, val{n});
			end
		end
		
		%	====================================================================

		function T = DefaultAnnealingScheduleFn(obj)
			T = 0.999*obj.Temperature;
		end
		
		%	====================================================================
		
		%	Randomly select co-ordintate of solution vector to change.
		%	The SolutionCost vector stores the cost value at the iteration in
		%	which each co-ordinate was changed. The probability of selection is
		%	proportional to the fractional cost for that index relative to all
		%	solution costs.
		function ind = DefaultIndexFn(obj)
			if length(obj.Soln)>2
				ind = randi(length(obj.Soln));
				while ind==obj.Index
					ind = randi(length(obj.Soln));
				end

			elseif length(obj.Soln)==2
				ind = 3 - obj.Index;
			else
				ind = 1;
			end
		end

		%	====================================================================

		function sol = DefaultNeighbourFn(obj)
			sol = obj.Solution;
			ind = obj.Index;
			sc = obj.StepSize(ind);
			sol(ind) = sol(ind).*(1 + sc) + sc;
		end
	
		%	====================================================================

		function val = DefaultAcceptanceFn(obj, CurCost)
			if CurCost==obj.Cost
				val = false;
			else
				val = CurCost<obj.Cost ...
					|| rand<exp((obj.Cost-CurCost)/obj.Temperature);
			end
		end
		
		%	====================================================================

		function val = DefaultConvergenceFn(obj)
			val = obj.Iteration>obj.IterMin ...
				&& (obj.Iteration>obj.IterMax ...
				|| obj.Temperature<=obj.MinTemperature);
		end
		
		%	====================================================================
		
		function DefaultPlotFn(obj)
			persistent h
			if obj.Iteration==obj.PlotRate || ~mod(obj.Iteration, 10*obj.PlotRate)
				N = 200;
				if obj.Iteration==obj.PlotRate
					h = [];
					x0 = 2*max(obj.Solution);
				else
					X = h.XData;
					Y = h.YData;
					x0 = 2*max(abs([X Y]));
				end
				x = linspace(-x0, x0, N);
				y = repmat(x, [N 1]);
				y = [Row(y); Row(y.')/2];
					imagesc(x, x/2, reshape(obj.CostFn(y), N, N));
				if isempty(h)
					h = line(obj.Solution(1), obj.Solution(2), ...
						'Marker', '+', 'Color', 'r');
				else
					h = line(X, Y, ...
						'Marker', '+', 'Color', 'r');
				end
				axis([-1 1 -.5 .5]*x0);
				grid on;
			end
			
			if obj.Iteration>obj.PlotRate
				x = h.XData;
				y = h.YData;
				if length(x)==100
					x = x(2:end);
					y = y(2:end);
				end
					
				h.XData = [x obj.Solution(1)];
				h.YData = [y obj.Solution(2)];
				mx = max([abs(x) abs(y)]);
				axis([-2 2 -1 1]*mx);

			end
			
			title(sprintf('%d: (x, y)_{min}=(%g, %g), (x, y)_{cur}=(%g, %g), T=%g, \\epsilon_{min}=%g,  \\epsilon_{cur}=%g', ...
				obj.Iteration, obj.OptimalSolution(1), obj.OptimalSolution(2), ...
				obj.Solution(1), obj.Solution(2), obj.Temperature, ...
				obj.OptimalCost, obj.Cost))
			drawnow;
		end
		
		%	====================================================================

	end %	methods (Access=private)
	
	
	
end %	classdef SimulatedAnnealing < handle
