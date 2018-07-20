function [varargout] = Swap(varargin)
%swap: Swaps multiple variables around
%
%	[y1, y2, ..., yn] = swap(x1, x2, ..., xn);
%
%	y1 = x1;
%	y2 = x2;
%	yn = xn;
%
%	Authors:
%		Dr Adam S Wyatt (a.wyatt1@physics.ox.ac.uk)
%
%	See also:

% for nv = 1:nargin
% 	varargout{nv} = varargin{nv};
% end

varargout = varargin(:);