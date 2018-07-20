function [dname, DATAFOLDER] = GetDataFolder(identifier, pth, AppendPath)
% Get data folder location based on current data analysis folder
%
%	[dname, DATAFOLDER] = GetDataFolder(identifier, pth, AppendPath)
%
%	dname		Data folder location
%	DATAFOLDER	Base data folder location
%	identifier	String to select which base folder to use [default = 20PWFE]
%	pth			(optional) current data analysis folder [default = pwd]
%	AppendPath	Append the data-analysis path to the base folder
%
%	This function will append all folders following 'data-analysis' in the
%	string pth (or current path) to the environment variable 'DATAFOLDER'
%

if ~exist('identifier', 'var') || isempty(identifier)
	identifier = '20PWFE';
end

if ~exist('pth', 'var') || isempty(pth)
	pth = pwd;
end

if ~exist('AppendPath', 'var') || isempty(AppendPath)
	AppendPath = true;
end

DATAFOLDERS = getenv('DATAFOLDER');
DATAFOLDER = regexpi(DATAFOLDERS, sprintf( ...
	'(?<=(%s[|])).*?(?=;)', identifier), 'match', 'once');

%	Add identifier if necessary
%	----------------------------------------------------------------------------
pth1 = fullfile(DATAFOLDER, identifier);
if exist(pth1, 'dir')
	DATAFOLDER = pth1;
end

%	Append input path
%	----------------------------------------------------------------------------
if AppendPath
	pth1 = regexpi(pth, '(?<=data-analysis).*', 'match', 'once');
	if isempty(pth1)
		dname = fullfile(DATAFOLDER, pth);
	else
		dname = fullfile(DATAFOLDER, pth1);
		if ~exist(dname, 'dir')
			dname = fullfile(DATAFOLDER, ...
				regexprep(pth1, identifier, '', 'once', 'ignorecase'));
% 			if ~exist(dname, 'dir')
% 				error('Could not find path');
% 			end
		end
	end
else
	dname = DATAFOLDER;
end

%	Check for duplicate strings
%	----------------------------------------------------------------------------
while ~exist(dname, 'dir')
	dname1 = regexprep(dname, '[\\/](\w+)[\\/]\1', [filesep '$1'], ...
		'once', 'ignorecase');
	if isequal(dname, dname1)
		error('Path not found');
	else
		dname = dname1;
	end
end

dname = regexprep(dname, '[\\/]', filesep);