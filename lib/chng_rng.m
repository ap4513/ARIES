function [y, c] = chng_rng(x, c)
%CHNG_RNG:	Converts between frequency and wavelength
%
%	y = chng_rng(x, c)
%
%	Uses formula:
%		y = 2*pi*c./x;
%
%	If c is not defined, c = 299.792458 is used (converts nm <--> rad/fs)
%
%  	Authors:
%  		Dr Adam S. Wyatt (a.wyatt1@physics.ox.ac.uk)

if (~exist('c', 'var') || isempty(c))
    c = 299.792458;
end

if iscell(x)
	y = cellfun(@(X) chng_rng(X, c), x, 'UniformOutput', false);
else
	y = CHNG(x, c);
end

end

function Y =  CHNG(X, c)
	Y = 2*pi*c./X;
end
