function [y0, var_m, w] = Weighted_Mean(y, w, Dim)
%Weighted_Mean: Performs a weighted mean
%
% Usage:
%
%	m	= Weighted_Mean(y, w, Dim)
%
%	m	= Weighted Mean
%	y	= Function to find mean
%	w	= Weights
%	Dim = Dimension to find mean over
%
%  	Authors:
%  		Dr Adam S. Wyatt (a.wyatt1@physics.ox.ac.uk)


%	Find 1st non-singular dimension
if (~exist('Dim', 'var')||isempty(Dim))
	Dim = find_first_dim(y);
end

Sy = sum(y.*w, Dim);

S = sum(w, Dim);
S(S==0) = eps;

y0 = Sy./S;

if nargout>=2
% 	var = sum((y-y0).^2.*w, Dim)/(Nx-1);
	var_m = sum((y-y0).^2.*w, Dim)./S;
	
	if nargout>=3
		w = w * (size(y, Dim)-1) ./ (S.*var_m);
	end
end