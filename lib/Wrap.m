function y = Wrap(x, interval, min)

if ~exist('min', 'var') || isempty(min)
	y = mod(x, interval);
else
	y = mod(x-min, interval) + min;
end

