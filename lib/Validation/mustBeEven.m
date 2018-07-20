function mustBeEven(val)

mustBeInteger(val);
if bitget(val, 1)
	error('Value must be even');
end

end
