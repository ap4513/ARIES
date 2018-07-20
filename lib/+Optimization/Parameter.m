classdef Parameter < handle & uiobjects.Object

	events
		ValueChanged;
		LowestChanged;
		HighestChanged;
	end
	
	%	========================================================================
	%	P R O P E R T I E S
	%	========================================================================
	
	properties
		SaveHistory logical;
		SaveHighest logical;
		SaveLowest logical;
		CoerceValues logical;
		NotifyChange logical;
	end %	properties
	
	%	========================================================================

	properties (SetAccess=protected)
		Value double;
		Initial double;
		Maximum double;
		Minimum double;
		Highest double;
		Lowest double;
	end
	
	%	========================================================================

	properties (Access=protected)
		HistoryLength double = 0;
		history double;
		UserData struct;
	end
	
	%	========================================================================

	properties (Dependent, SetAccess=public)
		History double;
	end	
	
	
	
	%	========================================================================
	%	M E T H O D S
	%	========================================================================
	
	methods

		%	====================================================================
		
		function obj = Parameter(varargin)
			
			%	Check if first argument is value (requires odd number of
			%	arguments)
			if bitget(nargin, 1)
				val = varargin{1};
				varargin(1) = [];
			else
				val = obj.FindArg(varargin, 'Value');
			end
			val = double(val);
			
			%	Set initial value (if not current value)
			obj.Initial = obj.FindArg(varargin, 'Initial', val);
							
			%	Set save flags
			obj.SaveHistory = obj.FindArg(varargin, 'SaveHistory', false);
			obj.SaveHighest = obj.FindArg(varargin, 'SaveHighest', false);
			obj.SaveLowest = obj.FindArg(varargin, 'SaveLowest', false);
			obj.NotifyChange = obj.FindArg(varargin, 'NotifyChange', false);

			%	Set minimum/maximum values
			obj.SetMaximum(obj.FindArg(varargin, 'Maximum', inf));
			obj.SetMinimum(obj.FindArg(varargin, 'Minimum', -inf));
			obj.CoerceValues = obj.FindArg(varargin, ...
				'CoerceValues', isfinite(obj.Maximum) || isfinite(obj.Minimum));

									
			%	Set lowest/highest values
			if obj.SaveHighest
				obj.Highest = obj.Initial;
			end
			
			if obj.SaveLowest
				obj.Lowest = obj.Initial;
			end
			
			%	Set user data
			UD = obj.FindArg(varargin, 'UserData');
			if ~isempty(UD)
				obj.UserData = struct(UD{:});
			else
				obj.UserData = struct;
			end
			
			%	Set value
			if ~isempty(val)
				obj.SetValue(val);
			else
				obj.SetValue(obj.Initial);
			end
		end %	function obj = Parameter(varargin)
		
		%	====================================================================
		
		function SetMaximum(obj, val)
			ind = (val<obj.Minimum);
			if any(ind)
				obj.Maximum(ind) = obj.Minimum(ind);
				obj.Minimum(ind) = val(ind);
			else
				obj.Maximum = val;
			end
		end %	function SetMaximum(obj, val)
		
		%	====================================================================
		
		function SetMinimum(obj, val)
			ind = (val>obj.Maximum);
			if any(ind)
				obj.Minimum(ind) = obj.Maximum(ind);
				obj.Maximum(ind) = val(ind);
			else
				obj.Minimum = val;
			end
		end %	function SetMinimum(obj, val)
		
		%	====================================================================
		
		function SetValue(obj, val)
			
			%	Coerce values
			%	----------------------------------------------------------------
			if obj.CoerceValues
				%	Set minimum
				ind = (val<obj.Minimum);
				if any(ind)
					if isscalar(obj.Minimum)
						val(ind) = obj.Minimum;
					else
						val(ind) = obj.Minimum(ind);
					end
				end
				
				%	Set maximum
				ind = (val>obj.Maximum);
				if any(ind)
					if isscalar(obj.Maximum)
						val(ind) = obj.Maximum;
					else
						val(ind) = obj.Maximum(ind);
					end
				end
			end
			
			%	Check if size of value has changed
			if ~isequal(size(obj.Value), size(val))
				obj.Reset;
			end
			
			%	Save extremeties
			%	----------------------------------------------------------------
			if obj.SaveHighest
				ind = (val>obj.Highest);
				obj.Highest(ind) = val(ind);
				if obj.NotifyChange
					notify(obj, 'HighestChanged');
				end
			end
			
			if obj.SaveLowest
				ind = (val<obj.Lowest);
				obj.Lowest(ind) = val(ind);
				if obj.NotifyChange
					notify(obj, 'LowestChanged');
				end
			end
			
			%	Save history
			%	----------------------------------------------------------------
			if obj.SaveHistory
				obj.AddToHistory(val);
			end
			
			%	Update value
			%	----------------------------------------------------------------
			obj.Value = val;
			if obj.NotifyChange
				notify(obj, 'ValueChanged');
			end
			
		end %	function SetValue(obj, val)
		
		%	====================================================================
		
		function val = get.History(obj)
			N = length(obj);
			if isvector(obj.Value)
				sz = [N obj.HistoryLength];
			else
				sz = [size(obj) obj.HistoryLength/N];
			end
			
			val = reshape(obj.history(1:obj.HistoryLength), sz);
		end %	function val = get.History(obj)
				
		%	====================================================================
		
		function Reset(obj)
			obj.HistoryLength = 0;
		end %	function Reset(obj)
		
		%	====================================================================
		
		function val = GetUserData(obj, str)
			if isfield(obj.UserData, str)
				val = obj.UserData.(str);
			else
				error('Field not found');
			end
		end
		
		%	====================================================================
		
		function SetUserData(obj, str, val)
			obj.UserData.(str) = val;
		end
		
		%	====================================================================
		
		function val = GetUserDataFieldNames(obj)
			val = fieldnames(obj.UserData);
		end
		
		%	====================================================================
		
		function val = double(obj)
			val = obj.Value;
		end
		
		%	--------------------------------------------------------------------
		
		function val = length(obj)
			val = numel(obj.Value);
		end
		
		%	--------------------------------------------------------------------

		function val = size(obj)
			val = size(obj.Value);
		end
		
		%	--------------------------------------------------------------------

		function val = ndims(obj)
			if isscalar(obj.Value)
				val = 0;
			elseif isvector(obj.Value)
				val = 1;
			else
				val = ndims(obj.Value);
			end
		end
		
		%	====================================================================
		
		function val = plus(obj1, obj2)
			val = double(obj1) + double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = minus(obj1, obj2)
			val = double(obj1) - double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = uminus(obj1)
			val = -double(obj1);
		end
		
		%	--------------------------------------------------------------------

		function val = uplus(obj1)
			val = double(obj1);
		end
		
		%	--------------------------------------------------------------------

		function val = times(obj1, obj2)
			val = double(obj1) .* double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = mtimes(obj1, obj2)
			val = double(obj1) * double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = rdivide(obj1, obj2)
			val = double(obj1) ./ double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = ldivide(obj1, obj2)
			val = double(obj1) .\ double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = mrdivide(obj1, obj2)
			val = double(obj1) / double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = mldivide(obj1, obj2)
			val = double(obj1) / double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = power(obj1, obj2)
			val = double(obj1) .^ double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = mpower(obj1, obj2)
			val = double(obj1) ^ double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = lt(obj1, obj2)
			val = double(obj1) < double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = gt(obj1, obj2)
			val = double(obj1) > double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = le(obj1, obj2)
			val = double(obj1) <= double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = ge(obj1, obj2)
			val = double(obj1) >= double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = ne(obj1, obj2)
			val = double(obj1) ~= double(obj2);
		end
		
		%	--------------------------------------------------------------------

		function val = colon(obj1, obj2, obj3)
			if nargin==2
				val = double(obj1):double(obj2);
			else
				val = double(obj1):double(obj2):double(obj3);
			end
		end
		
		%	--------------------------------------------------------------------

		function val = ctranspose(obj1)
			val = double(obj1)';
		end
		
		%	--------------------------------------------------------------------

		function val = transpose(obj1)
			val = double(obj1).';
		end
		
		%	--------------------------------------------------------------------

		function val = horzcat(varargin)
			val = cell2mat(cellfun(@(c) double(c), varargin, ...
				'UniformOutput', false));
		end
		
		%	--------------------------------------------------------------------

		function val = vertcat(varargin)
			val = cell2mat(cellfun(@(c) double(c), varargin.', ...
				'UniformOutput', false));
		end
		
		%	====================================================================

	end %	methods
	

	
	%	========================================================================
	%	M E T H O D S (Access=protected)
	%	========================================================================

	methods (Access=protected)
		
		%	====================================================================

		function GrowHistory(obj)
			if isempty(obj.history)
				obj.history = obj.Value;
			else
				obj.history = repmat(obj.history, [2 1]);
			end
		end %	function GrowHistory(obj)
		
		%	====================================================================

		function AddToHistory(obj, val)
			N = numel(val);
			ind = (1:N) + obj.HistoryLength;
			
			obj.HistoryLength = obj.HistoryLength + N;
			if obj.HistoryLength>length(obj.history)
				obj.GrowHistory;
			end
			
			obj.history(ind) = val;
			
		end %	function AddToHistory(obj, val)
	
		%	====================================================================

	end %	methods (Access=protected)
end %	classdef
