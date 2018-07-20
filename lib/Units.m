% ATOMIC_UNITS: Loads atomic unit definitions
%
%	USAGE:
%		[SI, AU] = Units;
%
%		value (SI) = value (AU) * AU.conversion
%
%
% Fundamental units
%	Atomic length (a0 = 5.3e-11m):
%		AU.a_0 = 5.29177210818e-11;
%	Mass of electron (me = 9.1e-31kg):
%		AU.m_e = 9.109382616e-31;
%	Electric charge (e = 1.602e-19C)
%		AU.e = 1.6021765314e-19;
%	Planck's constant (h/2pi = 1.05e-34Js)
%		AU.h_bar = 1.0545716818e-34;
%	Hartree energy (Eh = 4.36e-18J)
%		AU.En = 4.3597441775e-18;
%	Coulomb constant (1/4*pi*epsilon0 = 8.99e9Nm^2C^-2)
%		AU.Coulomb = 8.9875516e9;
%
% SI units
%	Boltzmann const (kB = 1.38e-23 J/K)
%		SI.k_B = 1.3806503e-23;
%	Permeability of free space (mu_0 = 1.26mkgs^-2A^-2)
%		SI.mu_0 = 1.25663706e-6;
%	Permittivity of free space (epsilon_0 = 8.85e-12s^4A^2m^-3kg^-1
%		SI.eps_0 = 8.85418782e-12;
%	Speed of light (c = 3e8m/s)
%		SI.C = 1/sqrt(SI.mu_0 * SI.eps_0);
%	Universal gas constant (J/K/mol):
%		SI.R = 8.314472;
%
% Other units
%	Time = 2.4e-17s
%		AU.time = AU.h_bar / AU.En;
%	Velocity = 2.2e6m/s
%		AU.velocity = AU.a_0 * AU.En / AU.h_bar;
%	Atomic force = 8.2e-8N
%		AU.force = AU.En / AU.a_0;
%	Current = 6.6e-3A
%		AU.current = AU.e * AU.En / AU.h_bar;
%	Boltzmann const in AU (kB_AU = 3.16e5K)
%		AU.k_B = AU.En / SI.k_B;
%	Pressure (P = 2.9e13Nm^-2)
%		AU.pressure = AU.En / AU.a_0^3;
%	Permeability of free space
%		AU.mu = AU.a_0 * AU.m_e / (AU.time^2 * AU.current^2);
%	Courtesy of Wiki: http://en.wikipedia.org/wiki/Atomic_units
%
%  	Authors:
%  		Dr Adam S. Wyatt (a.wyatt1@physics.ox.ac.uk)

function varargout = Units(varargin)
%% Fundamental units
%  -----------------------------------------------------------------------------

%	Atomic length (a0 = 5.3e-11m)
a_0 = 5.29177210818e-11;

%	Mass of electron (me = 9.1e-31kg)
m_e = 9.109382616e-31;

%	Electric charge (e = 1.602e-19C)
e = 1.6021765314e-19;

%	Planck's constant (h/2pi = 1.05e-34Js)
h_bar = 1.0545716818e-34;

%	Hartree energy (Eh = 4.36e-18J)
En = 4.3597441775e-18;

%	Coulomb constant (1/4*pi*epsilon0 = 8.99e9Nm^2C^-2)
Coulomb = 8.9875516e9;


%% SI units
%  -----------------------------------------------------------------------------

%	Boltzmann const (kB = 1.38e-23 J/K)
k_B_si = 1.3806503e-23;

%	Permeability of free space (mu_0 = 1.26mkgs^-2A^-2)
mu_0_si = 1.25663706e-6;

%	Permittivity of free space (epsilon_0 = 8.85e-12s^4A^2m^-3kg^-1
eps_0 = 8.85418782e-12;

%	Speed of light (c = 3e8m/s)
C = 1/sqrt(mu_0_si * eps_0);

%	Universal gas constant (J/K/mol):
R = 8.314472;

%% Other units
%  -----------------------------------------------------------------------------

%	Time = 2.4e-17s
time = h_bar / En;

%	Velocity = 2.2e6m/s
velocity = a_0 * En / h_bar;

%	Atomic force = 8.2e-8N
force = En / a_0;

%	Current = 6.6e-3A
current = e * En / h_bar;

%	Boltzmann const in AU (kB_AU = 3.16e5K)
k_B_au = En / k_B_si;

%	Pressure (P = 2.9e13Nm^-2)
pressure = En / a_0^3;

%	Permeability of free space
mu_0_au = a_0 * m_e / (time^2 * current^2);

%% Generate structures


if nargin==1
	if strcmpi(varargin, 'AU')
		varargout = {GetAU};
	elseif strcmpi(varargin, 'SI')
		varargout = {GetSI};
	else
		try
			varargout = {eval(varargin{1})};
		end
	end
elseif nargin
	varargout = cell(nargout, 1);
	for n=1:nargout
		try
			varargout{n} = eval(varargin{n});
		end
	end
else
	varargout = {GetSI, GetAU};
end

	function SI = GetSI
		SI = struct('k_B', k_B_si, 'mu_0', mu_0_si, 'eps_0', eps_0, ...
			'C', C, 'R', R);
	end
	
	function AU = GetAU
		AU = struct('a_0', a_0, 'm_e', m_e, 'e', e, 'h_bar', h_bar, ...
			'En', En, 'Coulomb', Coulomb, 'time', time, ...
			'velocity', velocity, 'force', force, 'current', current, ...
			'k_B', k_B_au, 'pressure', pressure, 'mu_0', mu_0_au);
	end
		
end