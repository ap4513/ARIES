%normalize Normalizes a vector/matrix to unit length/area
%
%	[A, A0] = normalize(A, [dim], [method], [x])
%
%	If dim==0, normalizes to whole data
%
%	switch(method)
%		case {2, 'peak'}
%			A0 = max(abs(A), [], dim);
%		case {1, 'intensity'}
%			A0 = sum(A, dim)
%		otherwise
%			A0 = sqrt(sum(abs(A).^2, dim))
%	end
%	A = A / A0

function [A, A0] = normalize(A, dim, method, x)
	
	if ~exist('dim', 'var') || isempty(dim)
		dim = find_first_dim(A);
		rsz = false;
	elseif dim==0
		dim = 1;
		rsz = true;
		N = size(A);
		A = A(:);
	else
		rsz = false;
	end
	
	if ~exist('method', 'var') || isempty(method)
		method = 0;
	elseif ischar(method)
			method = lower(method);
	end
	
	switch (method)
		case {2, 'peak'}
			A0 = max(abs(A), [], dim);
		case {1, 'intensity'}
			if exist('x', 'var') && ~isempty(x)
				A0 = abs(trapz(x, A, dim));
			else
				A0 = sum(A, dim);
			end
			
		otherwise
			if exist('x', 'var') && ~isempty(x)
				A0 = abs(realsqrt(trapz(x, abs(A).^2, dim)));
			else
				A0 = realsqrt(sum(abs(A).^2, dim));
			end
	end
	A = bsxfun(@rdivide, A, A0);

	if rsz
		A = reshape(A, N);
	end
	
end