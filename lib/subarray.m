function A = subarray(A, N, dims)
% SUBARRAY Extracts central subarray from larger array
%
%	A = subarray(A, N, [dims])
%
%	A		= Vector/array to select sub array from
%	N		= Vector of new dimension sizes
%	dims	= Vector of dimensions in which to extract from
%
%	Authors:
%		Dr Adam S Wyatt (a.wyatt1@physics.ox.ac.uk)
%
%	See also

if ~exist('dims', 'var') || isempty(dims)
	dims = (1:length(N));
end

sz = size(A);
[dims, ind] = unique(dims);
N = N(ind);

str = 'A = A(';
d = 1;
nd = length(dims);
for n=1:length(sz)
	if d<=nd && n==dims(d)
		M = floor((sz(dims(d))-N(d))/2);
		str = sprintf('%s%d:%d,', str, M+1, M+N(d));
		d = d+1;
	else
		str = [str ':,']; %#ok<*AGROW>
	end
end
str = [str(1:end-1) ');'];

eval(str);