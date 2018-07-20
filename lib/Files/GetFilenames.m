function fnames = GetFilenames(dname, str)
%	Obtain cellarray of filenames in folder matching a search string
%
%	fnames = GetFilenames(dname, str)
%
%	fnames	= Cell array of matched filenames
%	dname	= folder to search
%	str		= (optional) search string - default = all files
%
if ~exist('str', 'var') || isempty(str)
	fnames = dir(fullfile(dname, '*'));
else
	fnames = dir(fullfile(dname, str));
end

fnames = {fnames.name}.';
ind = ExtractStrings(fnames, '^\.{1,2}$');
fnames = fnames(~ind);
