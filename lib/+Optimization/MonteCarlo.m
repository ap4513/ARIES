classdef MonteCarlo < uiobjects.Object & handle
	
	properties
		Points double;
		Bounds double;
		Costs double;
		CostFun function_handle;
		Ndims double;
		Npoints double;
		MaxPoints double;
		TolX double = 1e-6;
		Iteration double = 0;
		NewPoints double;
		NewCosts double;
	end %	properties
	
	properties (Access=public)
		Bounds_ double;
		MinCost double;
		MaxCost double;
		MinPoint double;
		MaxPoint double;
	end
	
	methods
		
		function obj = MonteCarlo(fun, varargin)
			obj.CostFun = fun;
			
			bounds = obj.FindArg(varargin, 'Boundary');
			
			if isempty(bounds)
				bounds = 1;
			end
			
			if isscalar(bounds)
				
				N = obj.FindArg(varargin, 'Dimensions');
				if isempty(N)
					error('Must supply initial boundary or dimensionality');
				end
				bounds = ones(N, 1)*[-1 1]*bounds;
				
			elseif isvector(bounds)
				
				N = length(bounds);
				bounds = [-bounds bounds];

			else
				
				N = size(bounds);
				if length(N)>2 || N(2)~=2
					error('Bounds must be a vector [N, 1] or 2D matrix of size [N, 2]');
				end
				
			end
			
			obj.Ndims = N(1);
			obj.Bounds = sort(bounds, 2);
			obj.Iteration = 0;
			obj.Npoints = obj.FindArg(varargin, 'PointsPerIteration', 100);
			obj.MaxPoints = obj.FindArg(varargin, 'MaximumPoints', inf);
			
			obj.Bounds_ = [diff(bounds, [], 2) bounds(:, 1)];
		
		end %	functions
		
		function AddPoints(obj, N)
			if nargin<2
				N = obj.Npoints;
			end
			
			obj.Bounds_ = [diff(obj.Bounds, [], 2) obj.Bounds(:, 1)];
			points = obj.Bounds_(:, 2) ...
				+ (obj.Bounds_(:, 1) .* rand(obj.Ndims, N));
			
			if isempty(obj.Costs)
				obj.Costs = obj.CostFun(points);
				
				ind = ~isfinite(obj.Costs);
				obj.Costs(ind) = [];
				obj.Points = points(:, ~ind);
				
				[obj.MinCost, ind] = min(obj.Costs);
				obj.MinPoint = points(:, ind);
				
				[obj.MaxCost, ind] = max(obj.Costs);
				obj.MaxPoint = points(:, ind);
				
				obj.NewPoints = points;
				obj.NewCosts = obj.Costs;
			else
% 				obj.Points = [obj.Points points];
% 				obj.Costs = [obj.Costs obj.CostFun(points)];
% 				return

				%	Row: Scaled costs (Max..Min --> 0..1)
				weights = (obj.MaxCost - obj.Costs)/(obj.MaxCost - obj.MinCost).^2;
				weights = weights/sum(weights);

				%	Shift points to weighted centroid
				points = bsxfun(@plus, points, ...
					mean(bsxfun(@times, obj.Points, weights), 2));

				%	Scalar: RMS spread of points
				spread = mean(sum(abs(bsxfun(@minus, obj.Points, ...
					mean(obj.Points, 2))).^2, 1), 2);
				
				%	3D [x old new]: Vector new points --> old points
				Vec = bsxfun(@minus, obj.Points, reshape(points, [], 1, N));
				
				%	2D [1 old new]: Distance new points --> old points
				R2 = sum(abs(Vec).^2, 1);
				
				ind = squeeze(any(bsxfun(@le, R2, obj.TolX), 2));
				if any(ind)
					points(:, ind) = [];
					Vec(:, :, ind) = [];
					R2(:, :, ind) = [];
				end

				
				
				iweights = 1./weights;
				ind = ~isfinite(iweights);
				iweights(ind) = max(iweights(~ind(:)));
				iweights = iweights/sum(iweights);
				
				repulsion =  squeeze(sum(bsxfun(@times, Vec, ...
					bsxfun(@times, iweights, spread./(R2 + spread))), 2));
				attraction = squeeze(sum(bsxfun(@times, Vec, ...
					bsxfun(@rdivide, weights, R2)), 2));
				
				points = points + attraction - repulsion;
				ind = any(bsxfun(@gt, points, obj.Bounds(:, 2)), 1) ...
					| any(bsxfun(@lt, points, obj.Bounds(:, 1)), 1);
				points(:, ind) = [];
				costs = obj.CostFun(points);
				
				ind = ~isfinite(costs);
				if all(ind)
					return;
				end
				costs(ind) = [];
				points(:, ind) = [];
				
				[mx, ind] = max(costs);
				if mx>=obj.MaxCost
					obj.MaxCost = mx;
					obj.MaxPoint = points(:, ind);
				end
				
				[mn, ind] = min(costs);
				if mn<=obj.MinCost
					obj.MinCost = mn;
					obj.MinPoint = points(:, ind);
				end
				
				obj.NewPoints = points;
				obj.NewCosts = costs;
				
				obj.Points = [obj.Points points];
				obj.Costs = [obj.Costs costs];
				
				if size(obj.Points, 2)>obj.MaxPoints
					[costs, ind] = sort(obj.Costs);
					obj.Costs = costs(1:obj.MaxPoints);
					obj.Points = obj.Points(:, ind(1:obj.MaxPoints));
				end
				obj.Bounds = [max(obj.Points, [], 2)-min(obj.Points, [], 2) ...
					obj.Points(:, 1)];
					
			end

		end
		
	end %	methods
	
end %	classdef