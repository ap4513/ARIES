%	SINC
%		sinc(x) = sin(pi*x)/(pi*x)
function s = sinc(x)

ind = (x~=0);
s = ones(size(x));
x = pi*x(ind);
s(ind) = bsxfun(@rdivide, sin(x), x);