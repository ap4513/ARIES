%	Compute modulus squared
%
%	Usage:
%		r2 = abs2(z);
%
%	If z is real: r2 = z.^2
%	If z is complex: r2 = real(z).^2 + imag(z).^2
function r2 = abs2(z)

if isreal(z)
	r2 = z.^2;
else
	r2 = real(z).^2 + imag(z).^2;
end

end