function Et = Freq2Time(Ew, N, dim, shift)
%Freq2Time: Converts E(w) --> E(t)
%
%	Et = Freq2Time(Ew, [N], [dim], [shift])
%
%	Et		= Analytic field envelope
%	Ew		= Analytic spectrum
%	N		= Number of points (padded with zero if > size(Ew, dim)
%	dim		= Dimension to operate on
%	shift	= Perform (i)fftshift
%
%	shift can be a scalar or 2-element vector. If shift is a scalar, it is
%	converted to [shift; shift]. The default value is [0; 1]. If shift(1)~=0, an
%	ifftshift is performed prior to transforming. If shift(2)~=0, a fftshift is
%	performed after transforming. See NOTES below.
%
%	NOTES:
%	1)	(I)FFT assumes t(1) = 0, w(1) = 0.
%	2)	et = 2*real(Et.*exp(-i*w_dc*t)), where et is the real electric field
%		If min(w)>0, and spectral axis only just covers spectrum, real electric
%		field will be aliased - padding is required to resolve field.
%	3)	dw = (wmax-wmin)/(Nw-1) = 2*pi/(Nw*dt);
%		w = ((0:Nw-1)'-shift(1)*floor(Nw/2))*dw + wmin;
%		dt = 2*pi/(Nt*dw) = (tmax-tmin)/(Nt-1);
%		t = ((0:Nt-1)'-shift(2)*floor(Nt/2))*dt;
%	4)	if shift(1)==0, w_dc = w(1), else w_dc = w(floor(Nw/2)+1)
%		if shift(2)==0, t(1) = 0, else t(floor(Nt/2+1)) = 0
%	5)	smaller t = leading edge, longer t = trailing edge
%	6)	Dispersion: exp(i*phi(w)), normal dispersion ==> phi''>0 
%			==> red arrives at smaller t (blue larger t)
%
%	Authors:
%		Dr Adam S. Wyatt (a.wyatt1@physics.ox.ac.uk)
%
%	See also FFT, IFFT, TIME2FREQ, FREQ2TIMEQ, TIME2FREQQ


sz = size(Ew);

%	Check dim supplied
if (~exist('dim', 'var') || isempty(dim))
	dim = find_first_dim(Ew);
end

%	Check whether shift is applied
if (~exist('shift', 'var') || isempty(shift))
	shift = [0; 1];
elseif (numel(shift)==1)
	shift = [shift; shift];
end

%	Check Nt supplied
if (~exist('N', 'var') || isempty(N))
	N = sz(dim);
end


if N>sz(dim)
	Ew = padmatrix(Ew, N, dim);
elseif N<sz(dim)
	Ew = subarray(Ew, N, dim);
end

%	Perform shift (set w_dc index to floor(Nw/2)+1)
if shift(1)
	Ew = ifftshift(Ew, dim);
end

Et = fft(Ew, [], dim);

%	Perform shift (set t=0 index to floor(Nt/2)+1)
if shift(2)
	Et = fftshift(Et, dim);
end

end