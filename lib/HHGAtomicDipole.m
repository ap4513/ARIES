classdef HHGAtomicDipole < handle
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties
		%	Fundamental wavelength [nm]
		l1(1, 1) double {mustBePositive, mustBeFinite, mustBeNonempty} = 780;
		
		%	Harmonic orders
		q(:, 1) double {mustBePositive, mustBeFinite, mustBeNonempty} = 17;

		%	Nominal waist e^-2 radius [mm]
		%		This is the size of the waist (i.e. in the focus) of a Gaussian
		%		beam if there had not been any clipping
		w0(1,1) double {mustBePositive, mustBeFinite, mustBeNonempty} = .1;
		
		%	Iris location relative to nominal waist position [mm]
		zi(1,1) double {mustBeNegative, mustBeNonempty} = -750;
		%	Iris radius
		ai(1,1) double {mustBePositive, mustBeNonempty} = inf;
		%	Clipping resolution scale factor
		si(1,1) double {mustBePositive, mustBeNonempty, mustBeFinite} = 2;
		
		%	Target positions relative to nominal waist position [mm]
		z(:, 1) double {mustBeFinite, mustBeNonempty} = [-10:.2:0 0:.2:10].';
		
		%	Dipole phase co-efficient [10^-14 rad*cm^2/W]
		alpha_s(:,1) double {mustBeFinite, mustBeNonempty} = 4;
		alpha_l(:, 1) double {mustBeFinite, mustBeNonempty} = 22;

		%	Relative phase between short & long trajector [rad]
		%		dphi = phi0_s - phi0_l
		dphi(:, 1) double {mustBeFinite, mustBeGreaterThanOrEqual(dphi, ...
			-3.141592653589793115997963468544185161590576171875), ...
			mustBeLessThanOrEqual(dphi, ...
			3.141592653589793115997963468544185161590576171875)} = -pi/2;
		
		%	Effective nonlinearity
		%		I = |E1|^q_eff
		q_eff_s(:, 1) double {mustBeFinite, mustBePositive} = 7;
		q_eff_l(:, 1) double {mustBeFinite, mustBePositive} = 7;
		
		%	Trajectory intensity weighting: 
		%		I_s = beta * ...
		%		I_l = (1-beta) * ...
		beta(:, 1) double {mustBeFinite, mustBeNonnegative, ...
			mustBeLessThanOrEqual(beta, 1)} = .5;
		
		%	Distance from reference plane to detector plane [mm]
		L_det(1,1) double {mustBeFinite, mustBePositive} = 700;
		
		%	Peak intensity [10^14 W/cm^2]
		I0(1,1) double {mustBeFinite, mustBePositive} = 1;
	end

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (SetAccess=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess=private)

		%	Hankel matrix
		H = hankel_matrix(0, inf, 0);
		
		mu;
		
% 		IrrLq;	%	Interpolated XUV interference pattern
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (SetAccess=private, Dependent)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	%	Note: E(r) = U(r) .* obj.H.JR
	properties (SetAccess=private, Dependent)
		%	Dimension sizes
		Nq;		%	Number of harmonic orders
		Nz;		%	Number of target positions
		Nr;		%	Number of radial samples
		
		%	Harmonic wavelengths [nm]
		lq;		%	Harmonic wavelengths (Nq vec) [nm]
		
		%	Angular frequencies [rad/fs]
		w1;		%	Fundamental frequency (sc) [rad/fs]
		wq;		%	Harmonic frequencies (Nq vec) [rad/fs]
		
		%	Wave vector components [rad/mm]
		k1;		%	Fundamental wavenumber (sc) [rad/mm]
		kq;		%	Harmonic wavenumber (Nq vec) [rad/mm]
		kz1;	%	Fundamental longitudinal wave vector (Nr vec) [rad/mm]
		kzq;	%	Harmonic longitudinal wave vector (Nr x Nq) [rad/mm]

		%	Iris properties
		ri;		%	Radial co-ordinate at location of iris (sc) [mm]
		mi;		%	Aperture mask from iris (Nr vec)
		
		%	Fundamental fields
		Er0;	%	Nominal Gaussian field (Nr vec)
		Ur0;	%	Clipped Gaussian beam at nominal waist position (Nr vec)
		Urz1;	%	Clipped Gaussian beam around focus (Nr x Nz)
		
		%	Harmonic fields
		rL;		%	Radial co-ordinate at detector plane (Nr x Nq) [mm]
		Erzq;	%	Harmonic field around focus (Nr x Nz x Nq)
		Ur0q;	%	Harmonic field at nominal waist location (Nr x Nz x Nq)
		ErLq;	%	Harmonic intensity in far-field (Nr x Nz x Nq)
		IrLq;	%	XUV interference pattern (Nr x Nz/2 x Nz/2 x Nq)
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (Set/Get)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		%	Get number of harmonic orders
		function Nq = get.Nq(obj)
			Nq = length(obj.q);
		end
		
		%	Get number of target positions
		function Nz = get.Nz(obj)
			Nz = length(obj.z);
		end
		
		%	Get number of radial samples
		function Nr = get.Nr(obj)
			Nr = obj.H.Nr;
		end
		
		%	Harmonic wavelength [nm]
		function lq = get.lq(obj)
			lq = obj.l1./obj.q;
		end
		
		%	Fundamental frequency [rad/fs]
		function w1 = get.w1(obj)
			w1 = chng_rng(obj.l1);
		end
		
		%	Harmonic frequency [rad/fs]
		function wq = get.wq(obj)
			wq = obj.q .* obj.w1;
		end
		
		%	Fundamental wave number [rad/mm]
		function kl = get.k1(obj)
			kl = 2e6*pi./obj.l1;
		end
		
		%	Harmonic wave number [rad/mm]
		function kq = get.kq(obj)
			kq = obj.q .* obj.k1;
		end
		
		%	Fundamental longitudinal wave vector [rad/mm]
		function kz1 = get.kz1(obj)
			kz1 = obj.WaveVector(obj.k1);
		end
		
		%	Harmonic longitudinal wave vector [rad/mm]
		function kzq = get.kzq(obj)
			kzq = obj.WaveVector(Row(obj.kq));
		end
		
		%	Radial co-ordinate at iris [mm]
		function ri = get.ri(obj)
			ri = obj.FarFieldRadius(obj.k1, obj.zi);
		end
		
		%	Mask aperture
		function mi = get.mi(obj)
			delta = obj.si * diff(obj.ri(end-1:end));
			mi = .5*(1-tanh((obj.ri-obj.ai)/delta));
		end
		
		%	Gaussian field (normalized)
		function Er0 = get.Er0(obj)
			%	I(r) = exp(-2*(r/w0)^2) with w0 = e^-2 radius
			%	E(r) = sqrt(I(r)) (no phase for Gaussian at waist)
			Er0 = exp(-(obj.H.r./obj.w0).^2);
		end
		
		%	Clipped Gaussian beam
		function Ur0 = get.Ur0(obj)
			
			Ur0 = obj.Er0./obj.H.JR;
			
			%	Only propagate if iris has finite location and aperture
			if isfinite(obj.zi) && isfinite(obj.ai)
				phi = obj.FarFieldPhase(obj.k1, obj.zi);
				Uk = obj.FarField(Ur0, phi);
				Ur0 = obj.NearField(Uk .* obj.mi, conj(phi));
			end
		end

		%	Clipped beam through focus
		function Urz1 = get.Urz1(obj)
			Urz1 = obj.Propagate(obj.Ur0, obj.kz1, Row(obj.z));
		end
		
		%	Radial co-ordinate at detector [mm]
		function rL = get.rL(obj)
			rL = obj.FarFieldRadius(Row(obj.kq), obj.L_det);
		end
		
		%	Calculate harmonic field around focus
		function Erzq = get.Erzq(obj)
			Erzq = obj.HarmonicField(obj.Urz1.*obj.H.JR);
		end
		
		%	Calculate harmonic field at nominal waist location
		function Ur0q = get.Ur0q(obj)
			Ur0q = obj.Propagate(obj.Erzq./obj.H.JR, ...
				obj.Reshape(obj.kzq, obj.Nr), -Row(obj.z));
		end
		
		%	Calculate far-field harmonic field
		function ErLq = get.ErLq(obj)
			phi = obj.FarFieldPhase(obj.Reshape(obj.kq), obj.L_det);
			
			%	Note the far-field is effecitively in k-space, so should be
			%	scaled by J_v and not J_r
			ErLq = obj.FarField(obj.Ur0q, phi).*obj.H.JV;
		end
		
		%	Calculate far-field intensity pattern
		function IrLq = get.IrLq(obj)
			IrLq = obj.GetInterferogram(obj.ErLq);
		end
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		
		%	====================================================================
		%	Constructor
		%	====================================================================
		function obj = HHGAtomicDipole(varargin)
			%	Set up Hankel matrix
			%	----------------------------------------------------------------
			Nr = FindArg(varargin, 'RadialSamples', obj.Nr);%	Radial samples
			rmax = FindArg(varargin, 'MaxRadius', obj.H.R);	%	Max radius [mm]

			%	Check Hankel matrix is compatible with parameters
			if ~isequal(rmax, obj.H.R) || ~isequal(Nr, obj.Nr)
				%	Create Hankel matrix
				obj.H = hankel_matrix(0, rmax, Nr);
			end
			
			obj.si = FindArg(varargin, 'IrisResolutionScaling', obj.si);
			
			%	Set properties
			obj.SetProperties(varargin);
			
		end
		
		%	====================================================================
		%	Set object properties
		%	====================================================================
		%		These are the parameters that need to be minimized
		function SetProperties(obj, varargin)
			%	Central Wavelength [nm]
			obj.l1 = FindArg(varargin, 'Wavelength', obj.l1);

			%	Set harmonic order
			obj.q = FindArg(varargin, 'HarmonicOrders', obj.q);
			
			obj.w0 = FindArg(varargin, 'BeamWaist', obj.w0);

			%	Iris proprerties
			obj.zi = FindArg(varargin, 'IrisPosition', obj.zi);	
			obj.ai = FindArg(varargin, 'IrisRadius', obj.ai);

			%	Target positions
			obj.z = FindArg(varargin, 'TargetPositions', obj.z);
			
			%	Intensity dependent atomic dipole phase co-efficient
			obj.alpha_s = FindArg(varargin, 'AlphaShort', obj.alpha_s);
			obj.alpha_l = FindArg(varargin, 'AlphaLong', obj.alpha_l);
			
			%	Phase shift between short & long trajectory
			obj.dphi = FindArg(varargin, 'PhaseShift', obj.dphi);

			%	Effective nonlinearity
			obj.q_eff_s = FindArg(varargin, 'NonlinearityShort', obj.q_eff_s);
			obj.q_eff_l = FindArg(varargin, 'NonlinearityLong', obj.q_eff_l);
			
			%	Trajectory weighting
			obj.beta = FindArg(varargin, 'WeightingFactor', obj.beta);
	
			%	Detector distance
			obj.L_det = FindArg(varargin, 'DetectorDistance', obj.L_det);
			
			%	Peak intensity
			obj.I0 = FindArg(varargin, 'PeakIntensity', obj.I0);
		end
		
		%	====================================================================
		%	Get interpolated interferogram
		%	====================================================================
		function IrLq = InterpolatedInterferogram(obj, r, varargin)
			E = obj.ErLq;
			N = size(E);
			N(1) = length(r);
			Eq = zeros(N, 'like', 1i);
			rL = obj.rL;
			for nq=1:obj.Nq
				Eq(:, :, nq) = obj.interp(rL(:, nq), E(:, :, nq), r, ...
					varargin{:});
			end
			IrLq = obj.GetInterferogram(Eq);
		end
	
		%	====================================================================
		%	Calculate error between measured and simulated interferograms
		%	====================================================================
		function err = ErrorFun(obj, Im, r, varargin)
			
			RadialError = FindArg(varargin, 'RadialErrorDependence', true);
			RadialScaling = FindArg(varargin, 'RadialScaling', false);

			%	Remove background from measured data
			bg = FindArg(varargin, 'Background');
			if ~isempty(bg)
				Im = (Im-bg); %.*cast(Im>bg, 'like', Im);
			end

			%	Set XUV interference parameters
			obj.SetProperties(varargin);

			%	Simulate XUV interference patterns (interpolated onto MCP grid)
			method = FindArg(varargin, 'InterpMethod', 'makima');
			Is = obj.InterpolatedInterferogram(r, method{:});
% 			obj.IrrLq = Is;

			%	Calculate scaling factor (neglect saturated points)
			mu = FindArg(varargin, 'ScaleFactor');
% 			if ~isempty(mu)
% 				mu = reshape(mu, 1, 1, 1, []);
% 			end

			%	Calculate error
			%		If masked (saturated), then only include error if simulation
			%		is not saturated.
			mask = FindArg(varargin, 'Mask');
			
			%	Set scaling a function of (r,q) or (q)
			if RadialScaling
				SUM = @(x) sum(sum(x, 2), 3);
			else
				SUM = @(x) sum(sum(sum(x, 1), 2), 3);
			end
			
			if ~isempty(mask)
				nmask = ~mask;

				if isempty(mu)
					if RadialError
						mu = SUM(r.^2 .* Im.*Is.*mask) ...
							./ SUM(r.^2 .* Is.^2.*mask);
					else
						mu = SUM(Im.*Is.*mask)./SUM(Is.^2.*mask);
					end
				end
				Is = mu.*Is;
				
				%	Calcaulte differences
				if RadialError
					delta = r.*(Im - Is);
				else
					delta = (Im-Is);
				end
				
				%	Remove differences in which both measured and simulated
				%	values saturate
				ind = find(nmask);
				delta(ind(delta(nmask)>=0)) = [];
				
			else
				if isempty(mu)
					if RadialError
						mu = SUM(r.^2 .* Im.*Is)./SUM(r.^2 .* Is.^2);
					else
						mu = SUM(Im.*Is)./SUM(Is.^2);
					end
				end
				
				if RadialError
					delta = Column(r.*(Im - mu.*Is));
				else
					delta = Column(Im - mu.*Is);
				end
			end

			err = sqrt(mean(delta.^2));
			obj.mu = mu;
		end
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (Access=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access=private)
		
		%	====================================================================
		%	Calculate wave vector
		%	====================================================================
		function kz = WaveVector(obj, k)
			kz = k.^2 - obj.H.kr.^2;
			kz(kz<0) = 0;
			kz = sqrt(kz);
		end
		
		%	====================================================================
		%	Propagate field via ASM
		%	====================================================================
		%	Input and output fields are the scaled vectors, i.e. need to be
		%	divided/multiplied by obj.H.JR accordingly.
		function Urz = Propagate(obj, Ur0, kz, z)
			phi = exp(1i * kz .* z);
			Uk = obj.H.T * Ur0(:, :);
 			Urz = reshape(obj.H.T * ( Uk .* phi(:, :)), obj.Nr, length(z), []);
		end
		
		%	====================================================================
		%	Propagate field via far-field approximation
		%	====================================================================
		function r = FarFieldRadius(obj, k, z)
			r = abs(obj.H.kr .* (z ./ k));
		end
		
		function phi = FarFieldPhase(obj, k, z)
			phi = exp((.5i*k./z).*obj.H.r.^2);
		end
		
		function Urz = FarField(obj, Ur0, phi)
			U = Ur0.*phi;
			Urz = reshape(obj.H.T*U(:,:), ...
				obj.Nr, ...
				max(size(Ur0, 2), size(phi, 2)), ...
				max(size(Ur0, 3), size(phi, 3)));
		end
		
		function Ur0 = NearField(obj, Urz, phi)
			Ur0 = reshape((obj.H.T*Urz(:,:)).*phi, ...
				obj.Nr, ...
				max(size(Urz, 2), size(phi, 2)), ...
				max(size(Urz, 3), size(phi, 3)));
		end			
		
		%	====================================================================
		%	Calculate harmonic field
		%	====================================================================
		function Erzq = HarmonicField(obj, Erz1)
			
			%	Extract phase
			phiq = obj.Reshape(obj.q) .* angle(Erz1);
			
			%	Scale field
			I1 = obj.I0*abs2(Erz1);	%	Intensity
% 			I1 = obj.I0 * I1 ./ max(interp1(obj.H.r, I1, 0, 'makima'));
			A1 = sqrt(I1);		%	Amplitude

			%	Haromonic field
			%	----------------------------------------------------------------
			phis = phiq - obj.Reshape(obj.alpha_s).*I1 + obj.Reshape(obj.dphi);
			phil = phiq - obj.Reshape(obj.alpha_l).*I1;
			
			beta = obj.Reshape(obj.beta); %#ok<*PROPLC>
			Eqs = beta .* A1.^obj.Reshape(obj.q_eff_s) .* exp(1i*phis);
			Eql = (1-beta) .* A1.^obj.Reshape(obj.q_eff_l).* exp(1i*phil);

			Erzq = Eqs + Eql;
		end
		
		%	====================================================================
		%	Reshape function - sets harmonic order to 3rd dimension
		%	====================================================================
		function val = Reshape(~, x, sz1, sz2)
			
			%	By default, 1st & 2nd dimensions sizes are assumed to be unity
			%	Specify for non-unity dimension sizes
			if ~exist('sz1', 'var') || isempty(sz1)
				sz1 = 1;
			end
			if ~exist('sz2', 'var') || isempty(sz2)
				sz2 = 1;
			end
			
			%	Reshape to make harmonic order 3rd dimension
			val = reshape(x, sz1, sz2, []);
		end
		
		%	====================================================================
		%	Interpolation in complex plane
		%	====================================================================
		function E = interp(~, r, E, rr, varargin)
			rr = abs(rr);
			A = interp1(r, abs(E), rr, varargin{:});
			phi = interp1(r, unwrap(angle(E)), rr, varargin{:});
			E = A.*exp(1i*phi);
		end
		
		%	====================================================================
		%	Calculate interferogram
		%	====================================================================
		function Irzzq = GetInterferogram(~, Erzq)
			[Nr, Nz, Nq] = size(Erzq);
			N = Nz/2;
			Irzzq = abs2(reshape(Erzq(:, 1:N, :), Nr, N, 1, Nq) ...
				+ reshape(Erzq(:, N+1:end, :), Nr, 1, N, Nq));
		end
		
	end
	
end %	classdef
