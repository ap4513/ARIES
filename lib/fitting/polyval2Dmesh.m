function z = polyval2D(P, x, y)

Nx = length(x);
Ny = length(y);
[Qy, Qx] = size(P);

X = [Column(x).^(Qx-1:-1:1) ones(Nx, 1)].';
Y = [Column(y).^(Qy-1:-1:1) ones(Ny, 1)];

z = Y*P*X;
