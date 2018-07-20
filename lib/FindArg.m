% Find argument values from property name/value pairs
%
%	PropVal = FindARge(PropPairs, PropName, [DefaultValue])
%
%	PropVal:		Property value
%	PropPairs:		List of property value/pairs (e.g. varargin)
%	PropName:		Property name to find value of
%	DefaultValue:	(Optional) default value if property not found
%
%	Authors:	
%		Dr Adam S Wyatt (adam.wyatt@stfc.ac.uk)
%
%	See also FINDVAL
function PropVal = FindArg(PropPairs, PropName, DefaultValue)
	[~, PropVal] = FindVal(PropPairs, PropName);

	if isempty(PropVal) && exist('DefaultValue', 'var')
		PropVal = DefaultValue;
	end
end
