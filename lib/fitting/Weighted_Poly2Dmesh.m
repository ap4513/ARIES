function P = Weighted_Poly2D(x, y, z, wx, wy, Qx, Qy)

Nx = length(x);
Ny = length(y);

X = [x.^(Qx:-1:1) ones(Nx, 1)].';
Y = [y.^(Qy:-1:1) ones(Ny, 1)];

if exist('wx', 'var') && ~isempty(wx)
	wx = Row(wx);
	X = X.*wx;
	z = z.*wx;
end

if exist('wy', 'var') && ~isempty(wy)
	Y = Y.*wy;
	z = z.*wy;
end

P = Y\(z/X);

