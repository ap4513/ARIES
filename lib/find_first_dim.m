function [dim, sz] = find_first_dim(M)
%find_first_dim: Does exactly what it says on the tin
%
%	dim = find_first_dim(M)
%
%	Authors:
%		Dr Adam S Wyatt (a.wyatt1@physics.ox.ac.uk)

if numel(M)==1
	dim = 1;
	sz = 1;
else
	sz = size(M);
	dim = find(sz>1, 1, 'first');
	sz = sz(dim);
end
