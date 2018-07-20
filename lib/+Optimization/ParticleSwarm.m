classdef ParticleSwarm < uiobjects.Object & handle
	
	%	========================================================================
	%	P R O P E R T I E S
	%	========================================================================
	
	properties
		ErrorFun
		CurrentPosition
		
		DecayRate
		MutationRate
		CrossOverRate
		CrossOverDistance
		VelocityRate
		
		IterMin
		IterMax
				
		ErrorMin
		ErrorTol
		ConvergenceLength
		
		PlotOutput
		PlotRate
		PlotFun
		
		Repulsion;
		RepulsionOffset;
		
		Iter

	end %	properties
	
	%	========================================================================
	
	properties (SetAccess=private)
		Velocity
		Force

		CurrentError
		
		BestError
		BestPosition
		
		SolutionError
		Solution

		PosLength
		PosNum
		
		RelativeError
		
	end %	properties (SetAccess=private)
	
	%	========================================================================
	%	M E T H O D S
	%	========================================================================

	methods
		
		%	--------------------------------------------------------------------
		
		function obj = ParticleSwarm(InitPos, ErrorFun, varargin)
			
			obj.IterMin = obj.FindArg(varargin, 'IterMin', 0);
			obj.IterMax = obj.FindArg(varargin, 'IterMax', inf);
			obj.DecayRate = obj.FindArg(varargin, 'DecayRate', obj.IterMax/4);
			if ~isfinite(obj.DecayRate)
				obj.DecayRate = 100;
			end
			obj.MutationRate = obj.FindArg(varargin, 'MutationRate', .25);
			obj.CrossOverRate = obj.FindArg(varargin, 'CrossOverRate', .5);
			obj.ConvergenceLength = ...
				obj.FindArg(varargin, 'ConvergenceLength', 5);
			obj.ErrorTol = obj.FindArg(varargin, 'ErrorTol', 1e-6);
			obj.ErrorMin = obj.FindArg(varargin, 'ErrorMin', 1e-6);
			obj.PlotOutput = obj.FindArg(varargin, 'PlotOutput', false);
			obj.PlotRate = ...
				obj.FindArg(varargin, 'PlotRate', round(obj.IterMax/100));
			if ~isfinite(obj.PlotRate)
				obj.PlotRate = 1;
			end
			obj.PlotFun = obj.FindArg(varargin, 'PlotFun', @obj.Plot);
			obj.VelocityRate = obj.FindArg(varargin, 'VelocityRate', .75);
			obj.Repulsion = obj.FindArg(varargin, 'Repulsion', true);
			obj.RepulsionOffset = obj.FindArg(varargin, 'RepulsionOffset', .01);
			obj.CrossOverDistance = ...
				obj.FindArg(varargin, 'CrossOverDistance', .5);
			
			obj.CurrentPosition = InitPos;
			obj.ErrorFun = ErrorFun;
			[obj.PosLength, obj.PosNum] = size(obj.CurrentPosition);
			
			obj.Reset;
		end
		
		%	--------------------------------------------------------------------
		
		function Reset(obj)
			obj.CurrentError = obj.ErrorFun(obj.CurrentPosition);
			obj.BestPosition = obj.CurrentPosition;
			obj.BestError = obj.CurrentError;
			
			[~, ind] = min(obj.CurrentError);
			obj.Solution = obj.CurrentPosition(:, ind);
			obj.SolutionError = obj.CurrentError(ind);
			
			obj.RelativeError = inf(obj.ConvergenceLength, 1);
			
			obj.Iter = 0;
			obj.Velocity = zeros(obj.PosLength, obj.PosNum);
		end
		
		%	--------------------------------------------------------------------
	
		function Run(obj)
			
			Converged = false;
			
			while ~Converged
				
				T = (obj.Iter<obj.IterMin) + (obj.Iter>=obj.IterMin) ...
					*exp(-(obj.Iter-obj.IterMin)/obj.DecayRate);
				[mn, ind] = min(obj.CurrentError);
				if mn<obj.SolutionError
					obj.RelativeError = [obj.RelativeError(2:end); ...
						(obj.SolutionError-mn)/mn];

					obj.SolutionError = mn;
					obj.Solution = obj.CurrentPosition(:, ind);
				end
				
				ind = (obj.CurrentError<obj.BestError);
				obj.BestError(ind) = obj.CurrentError(ind);
				obj.BestPosition(:, ind) = obj.CurrentPosition(:, ind);
				
				%	Calculate force between current position and individual's
				%	best position
% 				dE1 = 1-exp(-abs((obj.CurrentError-obj.BestError) ...
% 					./obj.BestError))';
				dE1 = 1-exp(-abs(obj.CurrentError-obj.BestError)).';
				dP1 = obj.BestPosition-obj.CurrentPosition;
				N1 = (randn(size(obj.CurrentPosition))*obj.CrossOverRate*sqrt(T) ...
					+ obj.CrossOverDistance)*sqrt(T);
				F1 = bsxfun(@times, dE1, dP1).*N1;
				
				%	Calculate force between current position and global optimal
				%	position
% 				dE2 = 1 - exp(-abs((obj.CurrentError-obj.SolutionError) ...
% 					./obj.SolutionError))';
				dE2 = 1 - exp(-abs(obj.CurrentError-obj.SolutionError)).';
				dP2 = bsxfun(@minus, obj.Solution, obj.CurrentPosition);
				N2 = (randn(size(obj.CurrentPosition))*obj.CrossOverRate*sqrt(T) ...
					+ obj.CrossOverDistance);
				F2 = bsxfun(@times, dE2, dP2).*N2;
				
				%	Random mutation force
				if obj.Repulsion
					dP3 = obj.CurrentPosition./(obj.RepulsionOffset+dP2);
					dP3(~isfinite(dP3)) = 0;
				else
					dP3 = 1;
				end
				
				F3 = obj.MutationRate * T ...
					* randn(size(obj.CurrentPosition)).*dP3;
				
				obj.Force = F1 + F2 + F3;
				obj.Velocity = obj.VelocityRate*obj.Velocity + obj.Force;
% 				obj.Velocity = obj.Force;
				
% 				obj.CurrentPosition = obj.CurrentPosition + obj.Force;
				obj.CurrentPosition = obj.CurrentPosition + obj.Velocity;
				obj.CurrentError = obj.ErrorFun(obj.CurrentPosition);
				
				obj.Iter = obj.Iter + 1;
				
				if obj.PlotOutput && ~mod(obj.Iter, obj.PlotRate)
					obj.PlotFun(obj);
				end
				
				Converged = (obj.Iter>obj.IterMin && ...
					(obj.SolutionError<obj.ErrorMin ...
					|| all(obj.RelativeError<obj.ErrorTol) ...
					|| obj.Iter>obj.IterMax));

			end
			
		end
		
		%	--------------------------------------------------------------------

	end %	methods
	
	%	========================================================================
		
	methods (Access=private)
		function Plot(obj)
			X = linspace(-1, 1, obj.PosLength)';
			h = get(gca, 'Children');
			x = h(1).XData;
			ax = axis;

			plot(X, obj.CurrentPosition, ':', x, h(1).YData, '+-k');
			axis(ax);
			grid on;
			
			r = sort(obj.RelativeError);
			title(sprintf('Iteration: %d, e(%.3g, %.3g) = %.3g\ner = %s', ...
				obj.Iter, obj.Solution(1), obj.Solution(2), obj.SolutionError, ...
				num2str(r(1:5)', '%10.3g')));
			drawnow;
		end
	end %	methods (Access=private)
end
