%	Get the calling function filename and path
%		[pth, file] = GetFolder
%
%	Author(s):
%		Dr Adam S Wyatt (adam.wyatt@stfc.ac.uk)
function [pth, file] = GetFolder
	%	No input string - create default config file
	STFull = dbstack('-completenames');
	STFile = dbstack;

	pth = regexp(STFull(end).file, '^.*\\', 'match', 'once');
	file = STFile(end).file;
end
