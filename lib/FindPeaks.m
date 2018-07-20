function ind = FindPeaks(I, width, height, num)

%	Currently only works on vectors (which gets converted to a column)
if ~isvector(I)
	error('I must be a vector');
elseif isrow(I)
	I = I(:);
end

%	Check for width argument - default = 1
if ~exist('width', 'var') || isempty(width)
	width = 1;
end

N = length(I);

%	Array of "interior" points
ind0 = (width+1:N(1)-width).';

%	Array of "adjacent" points
ind1 = (-width:width);

%	Matrix of indices consisting of array of "adjacent" points for each
%	"interior" point
ind2 = ind0 + ind1;

%	Find all "interior" points greater than or equal to "adjacent" points
ind = find(all(I(ind0)>=I(ind2), 2)) + width;

%	If supplied, remove any values too low
if exist('height', 'var') && ~isempty(height)
	ind(I(ind)<height) = [];
end

%	If supplied, only return num largest points
if exist('num', 'var') && ~isempty(num) && num>0
	[~, indx] = sort(I(ind), 'descend');
	ind = ind(indx(1:min(num, length(ind))));
end
