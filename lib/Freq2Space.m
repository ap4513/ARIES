function Ex = Freq2Space(Ekx, N, dim, shift)
%Freq2Space: Converts E(w) --> E(t)
%
%	Ex = Freq2Space(Ekx, [N], [dim], [shift])
%
%	Needs editing
%	----------------------------------------------------------------------------
%
%	Ex		= Analytic field envelope
%	Ekx		= Analytic spectrum
%	N		= Number of points (padded with zero if > size(Ekx, dim)
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
%	2)	et = 2*real(Ex.*exp(-i*w_dc*t)), where et is the real electric field
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
%	----------------------------------------------------------------------------

sz = size(Ekx);

%	Check dim supplied
if (~exist('dim', 'var') || isempty(dim))
	dim = find_first_dim(Ekx);
end

%	Check whether shift is applied
if (~exist('shift', 'var') || isempty(shift))
	shift = [1; 1];
elseif (numel(shift)==1)
	shift = [shift; shift];
end

% 	Check Nt supplied
if (~exist('N', 'var') || isempty(N))
	N = sz(dim);
end

if (exist('N', 'var') && ~isempty(N))
	S = sz(dim);
	if N>S
		N = floor((N-S)/2);
		pad = zeros(1, length(sz));
		pad(dim) = N;
		Ekx = padarray(Ekx, pad);
	elseif N<S
		Ekx = subarray(Ekx, N, dim);
	end
end

%	Perform shift (set w_dc index to floor(Nw/2)+1)
if shift(1)
	Ekx = ifftshift(Ekx, dim);
end

Ex = ifft(Ekx, [], dim);

%	Perform shift (set t=0 index to floor(Nt/2)+1)
if shift(2)
	Ex = fftshift(Ex, dim);
end

end