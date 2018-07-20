%add_to_path: Adds folder & sub-folders to Matlab path, excluding .svn
%
%	Make sure your startup.m file points to this folder
%
%	Authors:
%		Dr Adam S Wyatt (a.wyatt1@physics.ox.ac.uk)
%
%	See also: startup
function output_path = add_to_path(current_dir)
	output_path = [];
	temp = regexp([genpath(current_dir) pathsep], ...
		['.[^' pathsep ']*' pathsep], 'match');
	for n=1:length(temp)
		if isempty(strfind(temp{n}, '.svn'))
			output_path = [output_path temp{n}];
		end
	end
end
