%% QUICK STRONG FIELD APPROXIMATION (QSFA)
% Simulates high harmonic generation (HHG) using a quick SFA model:
%	The action is calculated using the stationary momentum, but only 
%
%	Usage:
%		xt = QSFA(t, Et, Ip, tau)
%
%		xt = Time-dependent atomic dipole moment
%		t  = Temporal axis [au]
%		Et = Input real electric field [au]
%		Ip = Ionization potential [au]
%		tau_max = Maximum excursion time in continuum [au]
%		dipole = Dipole matrix function (Optional)
%		int_method = integration method (Optional)
%
%	dipole:
%		'@Hydrogenlike'		- default
%		'@Gaussian'
%		Handle to userdefined function, e.g. @MyFun
%
%	int_method:
%		'fft'
%		'sum'
%		'trapz'				- default
%%
function xt = QSFA(t, Et, Ip, tau_max, dipole, int_method)

if ~exist('dipole', 'var') || isempty(dipole)
	dipole = @Hydrogenlike;
elseif ischar(dipole)
	dipole = eval(dipole);
end

if ~exist('int_method', 'var') || isempty(int_method)
	int_method = 'trapz';
end

Nt = length(t);							%	Number of time samples
dt = diff(t(1:2));						%	Temporal resolution

DN = ceil(tau_max/dt);					%	Number of delay samples to integrate

%	Calculate indefinite integrals
switch lower(int_method)
	case 'fft'
		dw = 2*pi/(Nt*dt);
		w = ifftshift(((0:Nt-1).'-floor(Nt/2))*dw);
		iw = 1i./w;
		iw(~isfinite(iw)) = 0;

		Ew = Time2Freq(Et);
		Aw = -iw.*Ew;
		Bt = real(Freq2Time(iw.*Aw));
	% 	Bt = Bt - Bt(end);
		At = real(Freq2Time(Aw));
		At = At-At(end);
		B2t = real(.5*Freq2Time(Time2Freq(At.^2).*iw));

	case 'sum'
		At = -cumsum(Et)*dt;						%	Vector potential
		Bt = cumsum(At)*dt;						%	Inegral of VP
		B2t = .5*sum(At.^2)*dt;					%	Integral of .5*VP^2
		
	otherwise
		At = -cumtrapz(t, Et);					%	Vector potential
		Bt = cumtrapz(t, At);					%	Inegral of VP
		B2t = .5*cumtrapz(t, At.^2);			%	Integral of .5*VP^2
end

indr = (DN+1:Nt)';						%	Recombination indices
indi = (1:DN)' + (0:Nt-DN-1);			%	Ionization indices

tau = t(DN+1)-t(1:DN);					%	Delays

%	Calculate definite integrals
Bt = Bt(indr).' - Bt(indi);
B2t = B2t(indr).' - B2t(indi);


%	Stationary momentum
ps = Bt ./ tau;

%	Action
S = ((.5*ps.^2 + Ip) .* tau) + B2t - ps.*Bt;

%	Dipole moment at ionization time t_f
dpf = conj(dipole(ps - At(indr).', Ip));

%	Dipole moment at ionization time t_i
dpi = dipole(ps - At(indi), Ip);

A = (pi ./ (.5i*tau)).^1.5;				%	Amplitude factor

A = A .* (dpf .* dpi .* Et(indi) .* exp(-1i*S));

%	Dipole moment
xt = [zeros(DN, 1); dt*trapz(imag(A)).'];
end

function dp = Hydrogenlike(p, Ip)
alpha = 2*Ip;
A = 1i*2^3.5 * alpha^1.25 / pi;
dp = A * p .* (p.^2 + alpha).^-3;
end

function dp = Gaussian(p, Ip) %#ok<*DEFNU>
alpha = 2*Ip;
A = 1i*(pi*alpha).^-.75 ./ alpha;
dp = A * p .* exp(-p.^2/(2*alpha));
end