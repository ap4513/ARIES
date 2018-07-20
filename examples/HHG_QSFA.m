%% High Harmonic Generation via quick Strong Field Approximation

%	Scientific / Atomic units' constants
[SI, AU] = Units;
C = SI.C*1e3*1e-15;			%	Speed of light [mm/fs]

%	Pulse parameters
%	----------------------------------------------------------------------------
l0 = 780;									%	Central wavelength [nm]
w0 = chng_rng(l0);							%	Central frequency [rad/fs]
w0_eV = w0 * 1e15*AU.h_bar/AU.e;			%	Photon energy [eV]
T = 2*pi/w0;								%	Period [fs]
T_au = T * (1e-15/AU.time);					%	Period [au]
DT = 100;									%	Pulse duration, FWHM [fs]
DW = 4*log(2)/DT;							%	Pulse bandwidth, FWHM [rad/fs]
I0 = 4e14;									%	Peak intensity [W/cm^2]

%	Spectral/Temporal axes
%	----------------------------------------------------------------------------
Nw = 2^16;									%	Frequeny samples
qmax = 100;									%	Max harmonic order
wmax = qmax * w0;							%	Max frequency [rad/fs]
dw = 2*wmax/(Nw-1);							%	Frequency resolution [rad/fs]
dt = 2*pi/(Nw*dw);							%	Temporal resolution [fs]

w = (0:Nw-1)'*dw;							%	Frequency axis [rad/fs]
l = chng_rng(w);							%	Wavelength axis [nm]
t = ((0:Nw-1).'-floor(Nw/2))*dt;			%	Temporal axis [fs]
t_au = t*1e-15/AU.time;						%	Temporal axis [au]
w_eV = w * 1e15*AU.h_bar/AU.e;				%	Frequency axis [eV]
q = w/w0;									%	Harmonic order axis
k = w/C;									%	Wave vector [rad/mm]

%	Note the higher frequencies are simply the negative frequencies
wsq = (wmax-abs(w-wmax/2)).^2;		%	Freq axis squared

%	Fundamental pulse
%	----------------------------------------------------------------------------

%	Real electric field = 2*Re[E] where E is the analytic field - the factor of
%	2 is important!
Et = .5*gauss1D(t, 0, sqrt(2)*DT, 1).*exp(-1i*w0*t);
Ew = Time2Freq(Et);

%	Harmonic field
%	----------------------------------------------------------------------------

%	Ionisation potential
Ip_eV =	15.759;								%	Argon [eV]
Ip_au = Ip_eV * AU.e / AU.En;				%	[au]

%	Peak intensity and *real* field strength 
%		Note I = 1/2 eps_0 n c e^2 = 2 eps_0 n c |E|^2
%		Where e is the real field and E is the analytic field
e0 = sqrt(2e4*I0 / (SI.eps_0 * SI.C));		%	Peak real field strength [V/m]
e0_au = e0 / (AU.force / AU.e);				%	Peak real field strength [au]

%	Calcualte cut-off in eV
Up = 9.337*I0*1e-14*(l0*1e-3)^2;			%	Ponderomotive energy [eV]
wCO = Ip_eV + 3.17*Up;						%	Cut-off frequency [eV]

%	Harmonic spectrum at peak intensity
%		Note HHG strength is from acceleration of dipole strength, hence
%		requires multiplication of square of frequency
Ewq = Time2Freq(QSFA(t_au, e0_au*2*real(Et), Ip_au, T_au)) .* wsq;
Iwq = abs(Ewq).^2;

semilogy(q, Iwq);
line(wCO/w0_eV*[1 1], [min(Iwq)+eps max(Iwq)], ...
	'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
xlim([0 qmax]);
grid on