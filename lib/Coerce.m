%	Coerce values into defined range
%
%	[out, coerced] = Coerce(in, [mn], [mx]
%
%	If mn defined:
%		out(in<mn) = mn
%
%	If mx defined:
%		out(in>mx) = mx
%
%	Author(s)
%		Dr Adam S Wyatt (adam.wyatt@stfc.ac.uk)
function [out, coerced] = Coerce(in, mn, mx)

out = in;
% coerced = false(size(out));

if nargin>=2 && ~isempty(mn) && any(isfinite(mn(:)))
	out = max(out, mn);
end

if nargin>=3 && ~isempty(mx) && any(isfinite(mx(:)))
	out = min(out, mx);
end

if nargout>1
	coerced = (out~=in);
end

end
