function W = wigner(Ex)
%MYWIGNER: Calculates the Wigner distribution from a column vector
%
%	W  = wigner(Ex)
%
%	W  = output Wigner distribution
%	Ex = Input electric field (MUST be a column vector)
%
%	Notes:
%		W = Int(-inf..inf){E(x+y)E(x-y)exp[2ixy]}
%
%		E(x+y) & E(x-y) are calculated via a FFT (fast Fourier transform) using the
%		shift theorem. The integration is performed via a FFT. Thus it is important
%		for the data to satisfy the sampling theorem:
%		dy = 2*pi/X			X = span of all x-values	dy = y resolution
%		dx = 2*pi/Y			Y = span of all y-values	dx = x resolution
%		The data must be completely contained within the range x(0)..x(N-1) &
%		y(0)..y(N-1) (i.e. the function must fall to zero within this range).
%
%	v1.0
%
%	Currently waiting for update:
%		Allow an arbitrary output resolution
%		Allow an input vector for x (and possibly y).

%	Check data dimensions
N = size(Ex);
if any(N(2:end)-1)
    error('E must be a column vector');
end
N = N(1);

%	Generate shifting vectors
X = ((0:N-1)-N/2);
x = ifftshift(X'*2*pi/(N-1));
EXP = exp(bsxfun(@times, .5i*x, X));

%	Perform shift
Ek = fft(Ex);
EX1 = ifft(bsxfun(@times, Ek, EXP)); 
EX2 = ifft(bsxfun(@times, Ek, conj(EXP)));

W = real(fftshift(fft(ifftshift( ...
	bsxfun(@times, EX1.*conj(EX2), (abs(X)<=N/2)) , 2), [], 2), 2));