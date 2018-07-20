function P = polyfit2D(x, y, z, Qx, Qy)

Nx = length(x);
Ny = length(y);

X = [x.^(Qx:-1:1) ones(Nx, 1)].';
Y = [y.^(Qy:-1:1) ones(Ny, 1)];

P = Y\(z/X);

