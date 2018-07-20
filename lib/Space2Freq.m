function Ekx = Space2Freq(Ex, N, dim, shift)
%Space2Freq: Converts E(kx) --> E(x)
%
%	Ex = Space2Freq(Ekx, [N], [dim], [shift])
%
%	Ex		= Analytic field envelope
%	Ekx		= Analytic spectrum
%	N		= Number of points (padded with zero if > size(Ekx, dim)
%	dim		= Dimension to operate on
%	shift	= Perform (i)fftshift
%
%	shift can be a scalar or 2-element vector. If shift is a scalar, it is
%	converted to [shift; shift]. The default value is [1; 1]. If shift(1)~=0, an
%	ifftshift is performed prior to transforming. If shift(2)~=0, a fftshift is
%	performed after transforming. See NOTES below.
%
%	NOTES:
%	----------------------------------------------------------------------------
%	This needs editing
%	1)	(I)FFT assumes t(1) = 0, kx(1) = 0.
%	2)	ex = 2*real(Ex.*exp(-i*kx_dc*x)), where ex is the real electric field
%		If min(kx)>0, and spectral axis only just covers spectrum, real electric
%		field will be aliased - padding is required to resolve field.
%	3)	dkx = (kxmax-kxmin)/(Nkx-1) = 2*pi/(Nkx*dx);
%		kx = ((0:Nkx-1)'-shift(1)*floor(Nkx/2))*dkx + kxmin;
%		dx = 2*pi/(Nx*dkx) = (xmax-xmin)/(Nx-1);
%		t = ((0:Nt-1)'-shift(2)*floor(Nt/2))*dt;
%	4)	if shift(1)==0, w_dc = kx(1), else w_dc = kx(floor(Nw/2)+1)
%		if shift(2)==0, t(1) = 0, else t(floor(Nt/2+1)) = 0
%	5)	smaller t = leading edge, longer t = trailing edge
%	6)	Dispersion: exp(i*phi(kx)), normal dispersion ==> phi''>0 
%			==> red arrives at smaller t (blue larger t)
%
%	Authors:
%		Dr Adam S. Wyatt (a.wyatt1@physics.ox.ac.uk)
%
%	See also FFT, IFFT, FREQ2TIME, FREQ2TIMEQ, TIME2FREQQ
%	----------------------------------------------------------------------------


sz = size(Ex);

%	Check dim supplied
if (~exist('dim', 'var') || isempty(dim))
	dim = find_first_dim(Ex);
end

%	Check whether shift is applied
if (~exist('shift', 'var') || isempty(shift))
	shift = [1; 1];
elseif (numel(shift)==1)
	shift = [shift; shift];
end

%	Check Nt supplied
if (~exist('N', 'var') || isempty(N))
	N = sz(dim);
end

if (exist('N', 'var') && ~isempty(N))
	S = sz(dim);
	if N>S
		N = floor((N-S)/2);
		pad = zeros(1, length(sz));
		pad(dim) = N;
		Ex = padarray(Ex, pad);
	elseif N<S
		Ex = subarray(Ex, N, dim);
	end
end


%	Perform shift (set w_dc index to floor(Nw/2)+1)
if shift(1)
	Ex = ifftshift(Ex, dim);
end

Ekx = fft(Ex, N, dim);

%	Perform shift (set t=0 index to floor(Nt/2)+1)
if shift(2)
	Ekx = fftshift(Ekx, dim);
end

end