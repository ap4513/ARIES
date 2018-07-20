function [m, c, varm, varc, w] = Weighted_Linear(x, y, w, Dim)
% Weighted_Linear: Performs a weighted linear fit
%
%	[m, c, varm, varc, var] = Weighted_Linear(x, y, w, Dim)
%
%	m	= Slope
%	c	= Constant
%	x	= x-axis for y(x)
%	y	= Function to fit
%	w	= List of weights ~ sigma^-2
%	Dim = Dimension for x
%
%	min Sum[w*(y - (m*x + c))^2] for m, c
%
%  	Authors:
%  		Dr Adam S. Wyatt (a.wyatt1@physics.ox.ac.uk)


if (~exist('w', 'var')||isempty(w))
	w = ones(size(y));
end
if (~exist('Dim', 'var')||isempty(Dim))
	Dim = find(size(y)-1, 1, 'first');
end

S = sum(w, Dim);
Sx = sum(w.*x, Dim);
Sxx = sum(w.*x.^2, Dim);
Sy = sum(w.*y, Dim);
Sxy = sum(w.*x.*y, Dim);

%	D = (S.*Sxx - Sx.^2);
D = S.*Sxx - Sx.^2;

%	Calculate gradient
m = (S.*Sxy - Sx.*Sy)./D;

if nargout>=2
	%	Calculate offset
	c = (Sxx.*Sy - Sx.*Sxy)./D;
	
	if nargout>=3
		
		var = sum((y-(m.*x + c)).^2.*w, Dim)/(size(y, Dim)-2);
	
		%	Calculate variance in gradient
		varm = var.*S./D;
		
		if nargout>=4
			%	Calculate variance in offset
			varc = var.*Sxx./D;
			
			%	Calculate new weights
			if nargout>=5
				w = w./var;
			end
		end
	end
end

end
