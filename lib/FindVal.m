%	Find a list of property values from a list of property name/value pairs
%
%	[PropPairs, varargout] = FindProp(PropPairs, varargin)
%
%	PropPairs:		List of property name/value pairs with found pairs removed
%	varargout:		List of found property values
%	PropPairs:		List of property name/value pairs
%	varargin:		List of property names to search for
%
%	Example:
%		function out = MyFun(in, varargin)
%
%		[varargin, Val1, Val2] = FindProp(varargin, 'Name1', 'Name2');
%		if ~isempty(Val1)
%			...
%		end
%		...
%
%	Authors:	
%		Dr Adam S Wyatt (adam.wyatt@stfc.ac.uk)
%
%	See also FINDARG
function [PropPairs, varargout] = FindVal(PropPairs, varargin)
	if nargout<2
		varargout = {};
	else
		varargout = cell(nargout-1, 1);
	end

	while ~isempty(PropPairs) && length(PropPairs)<=1 && iscell(PropPairs{1})
		PropPairs = PropPairs{:};
	end

	%	Check for empty inputs
	if isempty(PropPairs) || isempty(varargin)
		return;
	end

	VarChck(PropPairs);

	% 	Loop over property names to search for
	for ni = 1:nargout-1
		%	Loop over property names to search from
		for np = 1:2:length(PropPairs)

			%	Save some typing
			p = PropPairs{np};
			v = varargin{ni};

			%	Check if start of property names are the same
			if strncmpi(v, p, min(length(v), length(p)))
				%	Extract property value
				varargout(ni) = PropPairs(np+1); %#ok<AGROW>

				%	Remove property pair from list
				PropPairs(np:np+1) = [];

				%	Go to next search for item
				break;
			end	%	if strncmpi( ...
		end %	for np = ...
	end %	for ni = ...
end	%	function ...


%	Check inputs
%		Ensures inputs form a set of property name/value pairs
function VarChck(var)
	pn = var(1:2:end-1);
	if mod(numel(var), 2)
		throwAsCaller(...
			MException( 'MATLAB:InvalidArgumentListLength', ...
			'Invalid argument list: must be PropName, PropVal pairs' ))
	elseif ~all( cellfun( @ischar, pn ) )
		throwAsCaller(...
			MException('MATLAB:InvalidArgument', ...
			'Invalid argument in list'))
	elseif any(arrayfun(@(n) any(strcmpi(pn(n), pn(n+1:end))), ...
			1:numel(pn)-1))
		throwAsCaller(MException('MATLAB:DuplicateArguments', ...
			'Duplicate property names in argument list'));
	end   
end % VarChk