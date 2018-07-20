%	Generate linear axes
%
%	Generates a linear uniformaly sampled frequency axis and its corresponding
%	(non-uniformly sampled) wavelength and (uniformly sampled) temporal axis.
%
%	Usage
%		[w, t1, dw, dt1, l, k, wmin, wmax, t, dt] = ...
%			GenerateAxes(Nw, lmin, lmax, {Nt1})
%
%	w		Linear uniformly sampled frequency axis [rad/fs]
%	t1		Linear uniformly sampled temporal axis [fs]
%	dw		Spectral resolution [rad/fs]
%	dt1		Temporal resolution [fs]
%	l		Wavelength axis [nm]
%	k		Wavevector axis [rad/nm]
%	wmax	Maximum frequency
%	wmin	Minimum frequency
%	t		Temporal axis (no padding: length(t1) = Nw) [fs]
%	dt		Temporal resolution (no padding: dt1=2*pi/(Nw*dw)) [fs]
%
%
%	Axes are suitable for use with Freq2Time & Time2Freq:
%
%		%	Initial parameters
% 		Nw = 2^10;
% 		lmin = 600;
% 		lmax = 1000;
% 		w0 = chng_rng(800);
% 		DT = 30;
%
%		%	Generate axes
% 		[w, t] = GenerateAxes(Nw, lmin, lmax);
%
%		%	Generate pulses
%		Iw = gauss1D(w, w0, 4*log(2)/DT, 1);
%		phiw = polyval([-3000, 300, 0, 0], w-w0);
%		Ew = sqrt(Iw).*exp(i*phiw); 
%		It0 = abs(Freq2Time(sqrt(Iw))).^2;
%		It = abs(Freq2Time(Ew)).^2;
%
%		%	Plot results
%		plot(t, abs(Freq2Time([abs(Ew) Ew])).^2); 
%		set(gca, 'XLim', [-300 300])
%
%	Author:
%		Dr Adam S Wyatt (adam.wyatt@stfc.ac.uk)

function [w, t1, dw, dt1, l, k, wmin, wmax, t, dt] = ...
	GenerateAxes(Nw, lmin, lmax, Nt)

wmax = chng_rng(lmin);
wmin = chng_rng(lmax);
w = linspace(wmin, wmax, Nw)';
dw = w(2)-w(1);

l = chng_rng(w);
k = 2*pi./l;

if ~exist('Nt', 'var') || isempty(Nt)
	Nt = Nw;
end

dt1 = 2*pi/(Nt*dw);
t1 = ((0:Nt-1)'-floor(Nt/2))*dt1;

dt = 2*pi/(Nw*dw);
t = ((0:Nw-1)'-floor(Nw/2))*dt;
