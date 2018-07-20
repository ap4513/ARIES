function mustBeOdd(val)

mustBeInteger(val);
if ~bitget(val, 1)
	error('Value must be odd');
end

end
