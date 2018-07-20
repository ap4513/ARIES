function M = padmatrix(M, Npad, dim)

%	Get matrix size
N = size(M);
N(dim) = Npad;
N = strrep(mat2str(N), ' ', ',');
str = ['M(' N(2:end-1) ') = 0;'];
eval(str);