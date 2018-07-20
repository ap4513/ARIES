classdef GeneticAlgorithm < handle
	
	%	========================================================================
	%	P R O P E R T I E S
	%	========================================================================
	
	properties
		ErrorFun
		ErrorMin
		ErrorTol

		Population
		PopNumMin
		PopAgeMax
		
		EliteNum

		DecayRate
		CrossOverRate
		MutationRate

		iterMax
		iterMin
		
		PopError
		OptimalError
		Solution
		
		ConvergenceLength
		PlotOutput
		
		PlotRate;

	end %	properties
	
	%	========================================================================
	
	properties (SetAccess=private)
		PopNum
		PopLength
		PopAge
		
		Elite
		Children
		
		iter
		RelError
		
	end %	properties (SetAccess=private)
	
	
	%	========================================================================
	%	M E T H O D S
	%	========================================================================
	
	methods
		
		%	--------------------------------------------------------------------
		
		function obj = GeneticAlgorithm(pop, fun, varargin)
			
			%	Get optional properties
			%	----------------------------------------------------------------
			obj.PopNumMin = FindArg(varargin, 'PopNumMin', 10);
			obj.iterMax = FindArg(varargin, 'IterMax', inf);
			obj.iterMin = FindArg(varargin, 'IterMin', 5);
			obj.ErrorMin = FindArg(varargin, 'ErrorMin', 0);
			obj.ErrorTol = FindArg(varargin, 'ErrorTol', 1e-6);
			obj.EliteNum = FindArg(varargin, 'EliteNum', 2);
			obj.PopAgeMax = FindArg(varargin, 'PopAgeMax', inf);
			obj.DecayRate = FindArg(varargin, 'DecayRate', obj.iterMax/4);
			obj.ConvergenceLength = FindArg(varargin, 'ConvergenceLength', ...
				max(obj.iterMin, obj.iterMax/10));
			obj.CrossOverRate = FindArg(varargin, 'CrossOverRate', .5);
			obj.MutationRate = FindArg(varargin, 'MutationRate', 1);
			obj.PlotOutput = FindArg(varargin, 'PlotOutput', false);
			obj.PlotRate = FindArg(varargin, 'PlotRate', 1);
			
			%	Initial conditions
			%	----------------------------------------------------------------
			obj.Population = pop;
			obj.ErrorFun = fun;
			
			[obj.PopLength, obj.PopNum] = size(obj.Population);

			%	Initialize iteration counter
			obj.Reset;
			
		end %	methods
		
		%	--------------------------------------------------------------------
		
		function Reset(obj)
			obj.PopAge = zeros(obj.PopNum, 1);
			obj.iter = 0;
			obj.OptimalError = inf;
			obj.Solution = nan(obj.PopLength, 1);
			obj.RelError = inf(obj.ConvergenceLength, 1);
			obj.PopError = ...
				[obj.ErrorFun(obj.Population(:, 1:obj.EliteNum)); ...
				zeros(obj.PopNum-obj.EliteNum, 1)];
		end
		
		%	--------------------------------------------------------------------
		
		function Run(obj)
			
			Converged = false;
			
			while ~Converged
				%	Get population error
				try
				obj.PopError = [obj.PopError(1:obj.EliteNum); ...
					obj.ErrorFun(obj.Population(:, obj.EliteNum+1:end))];
				catch
					fprintf('Oops\n');
				end
				
				%	Sort errors
				[obj.PopError, ind] = sort(obj.PopError);

				%	Get optimal solution
				if obj.PopError(1)<obj.OptimalError
					obj.RelError = [obj.RelError(2:end); ...
						(obj.OptimalError-obj.PopError(1)) ...
						/obj.PopError(1)];
					obj.OptimalError = obj.PopError(1);
					obj.Solution = obj.Population(:, ind(1));
				end
				
				%	Extract elite populations
				obj.Elite = obj.Population(:, ind(1:obj.EliteNum));
				
				%	Kill-off 'oldies'
				obj.Population(:, (obj.PopAge>obj.PopAgeMax)) = [];
				obj.PopNum = size(obj.Population, 2);

				%	Crossover and mutation
				obj.CrossOver;
				obj.Mutation;
				
				%	Generate new population
				obj.Population = [obj.Elite obj.Children];
				obj.PopAge = [obj.PopAge(ind(1:obj.EliteNum))+1; ...
					zeros(size(obj.Children, 2), 1)];
				obj.PopNum = size(obj.Population, 2);
				
	
% 				obj.Population = unwrap(obj.Population);
% 				obj.Population = bsxfun(@minus, obj.Population, mean(obj.Population));
				obj.iter = obj.iter + 1;
				
				Converged = (obj.iter>obj.iterMin && ...
					(obj.OptimalError<obj.ErrorMin || ...
					all(obj.RelError<obj.ErrorTol) || obj.iter>obj.iterMax));
				
				if obj.PlotOutput && ~mod(obj.iter, obj.PlotRate)
					obj.Plot;
				end
			end
			
		end
		
		%	--------------------------------------------------------------------
		
	end
	
	%	========================================================================
	
	methods (Access=private)
		
		function CrossOver(obj)
			N = max(obj.PopNum, obj.PopNumMin) - obj.EliteNum;
			
			%	Select children (binary tournement)
			ind1 = randi([1 obj.PopNum], N, 1);
			ind2 = randi([1 obj.PopNum], N, 1);
			ind3 = obj.PopError(ind1)<obj.PopError(ind2);
			obj.Children = obj.Population(:, [ind1(ind3); ind2(~ind3)]);
			
			ind1 = randi([1 N], N, 1);
			ind2 = randi([1 N], N, 1);
			Delta = obj.Children(:, ind1) - obj.Children(:, ind2);
			obj.Children = obj.Children ...
				+ randn(obj.PopLength, N).*Delta*obj.CrossOverRate;

% 			%	Select child pairs
% 			M = floor(N/2);
% 			ind1 = randi([1 N], M, 1);
% 			ind2 = randi([1 N], M, 1);
% 
% 			%	Random numbers
% % 			ind = (rand(1, M)>.5);
% 			beta = obj.CrossOverRate*log(rand(1, M));
% % 			beta(ind) = -beta(ind);
% 			
% 			s = bsxfun(@times, beta, ...
% 				abs(obj.Children(:, ind1) - obj.Children(:, ind2)));
% 			obj.Children = obj.Children(:, [ind1; ind2]) + [s s];
% 			C1 = obj.Children(:, ind1) 
% 			C2 = obj.Children(:, ind2) ...
% 				+ beta.*Delta;	
% 			obj.Children = [C1 C2];
		end
		
		%	--------------------------------------------------------------------
		
		function Mutation(obj)
			%	'Temperature' scaling factor
			if obj.iter>obj.DecayRate
				return;
			end
			
			N = size(obj.Children, 2);
			T = (1 - rand(obj.PopLength, N).^(1 - obj.iter/obj.DecayRate))...
				.^obj.MutationRate;
			

			%	Population selection
			ind = (rand(obj.PopLength, N)<=.5);
			
			%	Mutate 1st half
			Delta = bsxfun(@minus, max(obj.Children), obj.Children);
			obj.Children(ind) = obj.Children(ind) ...
				+ Delta(ind).*T(ind);

			%	Mutate 2nd hald
			ind = (~ind);
			Delta = bsxfun(@minus, obj.Children, min(obj.Children));
			obj.Children(ind) = obj.Children(ind) ...
				- Delta(ind).*T(ind);
			
		end
		
		%	--------------------------------------------------------------------

		function Plot(obj)
			ax = axis;
			x = linspace(-1, 1, obj.PopLength)';
			h = get(gca, 'Children');
			plot(x, (obj.Population), ':', ...
				h(1).XData, [spline(x, obj.Solution, h(1).XData(:)) h(1).YData(:)], '+-');
% 			axis tight
% 			axis([-1 1 -2 2])
			axis(ax);
			title([obj.iter obj.PopNum toc obj.OptimalError max(obj.RelError)]);
			drawnow;
		end
		
		%	--------------------------------------------------------------------
		
	end %	methods (Access=private)
	
	%	========================================================================
	
end
