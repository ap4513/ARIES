%arrow2: Draws a 3D arrow
%
%	h = arrow2;
%	h = arrow2(arrow2);
%	h = arrow2(start, finish, radii, lengths, color, edges);
%	h = arrow2(..., Property1, Value1, Property2, Value2, ...);
%	h = arrow2(Property1, Value1, Property2, Value2, ...);
%
%	radii:
%		Set the first and last values to zero to close ends.
%		For standard arrow, use radii = [0 r1 r2 r2 r2 0], where r1 and r2 are the
%		tail and head radii respectively
%
%	lengths:
%		The first and last elements are always 0 and 1 respectively.
%		The relative position of each radii
%
%	color:
%		Can either be a flat color or CData
%
%	edges:
%		true/false - determines whether to draw edges
%
%	Author:
%		Adam S Wyatt

classdef arrow2 < handle

	%	========================================================================
	
	properties
		start;
		finish;
		radii;
		lengths;
		color;
		handle;
		edges;
	end
	
	%	========================================================================

	properties (SetAccess = protected)
		theta;
		X;
		Y;
		Nr;
	end
			
	%	========================================================================

	methods

		%	--------------------------------------------------------------------
		
		%	Constructor
		function obj = arrow2(varargin)
			Init = {'start', [0 0]; ...
				'finish', [1 0]; ...
				'radii', [0 .1 .1 .5 0]; ...
				'lengths', [0 0 .75 .75 1]; ...
				'color', [0 0 1]; ...
				'edges', false};
			
			%	Check input type
			if nargin<1
				
				%	Empty - use defaults
				for n=1:size(Init, 1)
					obj.(Init{n, 1}) = Init{n, 2};
				end
				
				Properties = {};
				Values = {};
				
			elseif isa(varargin{1}, 'arrow2')

				Properties = varargin(2:2:end-1);
				Values = varargin(3:2:end);
				
				old = varargin{1};
				for n=1:size(Init, 1)
					obj.(Init{n, 1}) = old.(Init{n,1});
				end
				
			elseif ischar(varargin{1})
				
				%	Property/Value pair
				Properties = varargin(1:2:end-1);
				Values = varargin(2:2:end);
				for n=1:size(Init, 1)
					ind = strcmpi(Properties, Init{n,1});
					if any(ind)
						obj.(Init{n,1}) = Values{ind};
						Properties(ind) = [];
						Values(ind) = [];
					else
						obj.(Init{n,1}) = Init{n,2};
					end
				end
				
			else

				%	Ordered list
				for N=1:nargin
					if ischar(varargin{N})
						N = N-1;
						break;
					end
				end
				
				for n=1:N
					if ~isempty(varargin{n})
						obj.(Init{n,1}) = varargin{n};
					else
						obj.(Init{n,1}) = Init{n,2};
					end
				end
				
				for n=N+1:size(Init, 1)
					obj.(Init{n, 1}) = Init{n, 2};
				end
				
				for n=nargin+1:size(Init,1)
					obj.(Init{n,1}) = Init{n, 2};
				end
				
				Properties = varargin(N+1:2:end-1);
				Values = varargin(N+2:2:end);
				
			end %	Input type
			
			obj.setcoord;
			obj.draw;
			for n=1:length(Properties)
				set(obj.handle, Properties{n}, Values{n});
			end
			
		end %	Constructor
	
		%	--------------------------------------------------------------------
		
		function set.start(obj, s)
			if numel(s)~=2
				error('Start must be [x, y] vector')

			else
				obj.start = s(:);
			end
			
			if ~isempty(obj.X)
				obj.setcoord;
			end
		end

		%	--------------------------------------------------------------------
		
		function set.finish(obj, f)
			if numel(f)~=2
				error('Finish must be [x, y] vector');
			else
				obj.finish = f(:);
			end
			
			if ~isempty(obj.X)
				obj.setcoord;
			end
		end
		
		%	--------------------------------------------------------------------
		
		function set.radii(obj, r)
			
			if any(r<0)
				error('Radii must be greater than or equal to zero.');
			end
			
			if length(r)<2
				error('Radii must contain at least 2 elements.');
			end
			
			obj.Nr = numel(r);
			obj.radii = reshape(r, 1, []);
			
			if ~isempty(obj.X)
				obj.setcoord;
			end
		end
				
		%	--------------------------------------------------------------------
		
		function set.lengths(obj, l)
			Nl = length(l);
			
			if any(l<0) || any(l>1)
				error('Lengths must be in range (0,1)');
			end
			
			if Nl~=obj.Nr
				if Nl~=(obj.Nr-2)
					error('length(lengths) must be equal to length(radii)');
				else
					l = [0; l(:); 1];
				end
			else
				if l(1) ~= 0
					l(1) = 0;
				end
				if l(end) ~= 1;
					l(end) = 1;
				end
			end
			
			obj.lengths = reshape(l, 1, []);

			if ~isempty(obj.X)
				obj.setcoord;
			end
		end
				
		%	--------------------------------------------------------------------

		function set.color(obj, c)
			if size(c)==[1, 3]
				obj.color = c;
				if ~isempty(obj.handle)
					set(obj.handle, 'FaceColor', c);
				end
			elseif size(c)==size(obj.Z)
				obj.color = c;
				if ~isempty(obj.handle)
					set(obj.handle, 'CData', c);
				end
			end
		end
			
	end %	methods
	
	%	========================================================================

	methods (Access = private)

		%	--------------------------------------------------------------------
		
		function setcoord(obj)
			delta = obj.finish-obj.start;
			L = realsqrt(sum(delta.^2));
			theta = atan2(delta(2), delta(1));
			
			V = [L*[obj.lengths fliplr(obj.lengths(1:end-1))]; ...
				[obj.radii -fliplr(obj.radii(1:end-1))]];
			V = arrow2.Rxy(theta)*V;
			
			obj.X = obj.start(1) + V(1, :).';
			obj.Y = obj.start(2) + V(2, :).';
			
			if ~isempty(obj.handle)
				obj.redraw;
			end
		end

		%	--------------------------------------------------------------------

		function redraw(obj)
			set(obj.handle, 'XData', obj.X, 'YData', obj.Y);
			if size(obj.Y)==obj.color
				set(obj.handle, 'CData', obj.color);
			else
				set(obj.handle, 'FaceColor', obj.color);
			end
		end

		%	--------------------------------------------------------------------

		function draw(obj)
			obj.handle = patch(obj.X, obj.Y, obj.color);

			if size(obj.color)==[1, 3]
				set(obj.handle, 'FaceColor', obj.color);
			end
			
			if ~obj.edges
				set(obj.handle, 'Linestyle', 'none');
			end
			
		end

		%	--------------------------------------------------------------------
	
	end
	
	%	========================================================================

	methods (Static)
		
		%	--------------------------------------------------------------------
		
		function r = Rxy(theta)
			r = [	cos(theta) -sin(theta); ...
					sin(theta)	cos(theta)];
		end
		
		%	--------------------------------------------------------------------
		
	end %	methods (Static)
	
	%	========================================================================

end %	classdef