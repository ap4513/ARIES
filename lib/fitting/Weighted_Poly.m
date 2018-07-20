%% Weighted_Poly: Performs weighted polynomial fit
%
%	P	= Weighted_Poly(x, y, w, n)
%
%	P	= Column vector of polynomial co-efficients (n:-1:0)'
%
%	x	= x-axis (column vector)
%	y	= Data (matrix of column vectors)
%	w	= Weights (e.g. 1./sigma(y) - column vector)
%	n	= Order of fitting
%
%	Uses ldivide to calculate the weighted fit of the function:
%		Y = sum(P(i+1) * xi^(n-i), i=0..n)
%
%  	Authors:
%  		Dr Adam S. Wyatt (a.wyatt1@physics.ox.ac.uk)
%
%	See also LDIVIDE
%
%	Numerical Receipes 3rd Edition (15.4.2)
function P = Weighted_Poly(x, y, w, n)	

A = bsxfun(@power, x, n:-1:0);
if exist('w', 'var') && ~isempty(w)
	A = bsxfun(@times, A, w);
	y = bsxfun(@times, y, w);
end

P = A\y;

