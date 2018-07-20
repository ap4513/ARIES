classdef RasterScan < handle
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties
		RasterFolder(1, :) char = '2017-08-22_190845_Sort_focus_345uJ_raster';
		
		%	MCP spatial resolution [mm]
		dr_MCP(1,1) double {mustBePositive, mustBeFinite} = 44.48e-3;
		
		
		%	On-axis pixel location
		nth0(1,1) double {mustBePositive, mustBeFinite} = 1;

		%	On-harmonic pixel locations
		nw0(:, 1) double {mustBePositive, mustBeFinite, mustBeNonempty} = 1;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (SetAccess=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess=private)
		BaseFolder(1, :) char = 'C:\Data';
		ScanData(:, 1) cell = {};
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (SetAccess=private, Dependent)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess=private, Dependent)
		Folder;
		CamInfo;
		ScanInfo;
		
		zT1;			%	Target 1 stage positions [mm]
		zT2;			%	Target 2 stage positions [mm]
		NzT1;			%	Number of T1 steps
		NzT2;			%	Number of T2 steps

		Nq;				%	Number of harmonic orders
		Nth;			%	Height of image (divergence) [pxl]
		Nw;				%	Width of image (frequency) [pxl]
		nth;			%	Divergence axis [pxl]
		nw;				%	Frequency axis [pxl]

	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		
		%	====================================================================
		%		C O N S T R U C T O R
		%	====================================================================
		function obj = RasterScan(varargin)
			%	Set base folder
			DefaultBaseFolder = {'C:\Data'; 'D:\Data'; ...
				'\\physics\dfs\DAQ\Atomic & Laser\Hooker'};
			DefaultSubFolder = { ...
				'Artemis\Experiment data\17120006 OKeeffe\Raw Data'; ...
				'Artemis 2017\17120006 OKeeffe\Raw Data'};
			obj.SetBaseFolder( ...
				FindArg(varargin, 'BaseFolder', DefaultBaseFolder), ...
				FindArg(varargin, 'SubFolder', DefaultSubFolder));
			
			%	Set raster scan folder
			obj.RasterFolder = FindArg(varargin, 'RasterScan', ...
				obj.RasterFolder);
			
			%	Load raster scan data
			LoadRasterScan(obj);
			
			%	Set parameters
			obj.dr_MCP = FindArg(varargin, 'MCPResolution', obj.dr_MCP);

		end
		
		%	====================================================================
		%		S E T   B A S E   F O L D E R
		%	====================================================================
		function SetBaseFolder(obj, BaseDir, SubBaseDir)
			%	Find which of these folders actually exists
			if iscell(BaseDir)
				n = find(cellfun(@(str) exist(str, 'dir'), BaseDir), 1, 'first');
				if isempty(n)
					error('Could not find data location!');
				end
				BaseDir = BaseDir{n};
			end
			
			if nargin>2 && ~isempty(SubBaseDir)
				if iscell(SubBaseDir)
					%	Find which of these folders actually exists
					n = find(cellfun(@(str) exist(fullfile(BaseDir, str), ...
						'dir'), SubBaseDir), 1, 'first');
					if isempty(n)
						error('Could not find data sub folder location!');
					end

					%	Combine base folder, sub-folder and data folder
					BaseDir = fullfile(BaseDir, SubBaseDir{n});
				else
					BaseDir = fullfile(BaseDir, SubBaseDir);
				end
			end
			
			obj.BaseFolder = BaseDir;
		end
		
		%	====================================================================
		%		L O A D   R A S T E R   S C A N
		%	====================================================================
		function obj = LoadRasterScan(obj, folder)
			%	Check if optional raster scan folder is supplied
			if nargin>1 && ~isempty(folder)
				obj.RasterFolder = folder;
			end
			
			%	Load scan data
			fnames = ExtractStrings(GetFilenames(obj.Folder), 'ScanData', true);
			obj.ScanData = cellfun(@(str) ...
				matfile(fullfile(obj.Folder, str), 'Writable', false), ...
				fnames, 'UniformOutput', false);
		end
		
		%	====================================================================
		%		G E T   D A T A
		%	====================================================================
		function I_r_z1_z2_w = GetData(obj, indth, indz1, indz2, indw, fast)
			%	Check inputs
			%	----------------------------------------------------------------
			if nargin<2 || isempty(indth)
				indth = obj.nth;
			end
			if nargin<3 || isempty(indz1)
				indz1 = 1:obj.NzT1;
			end
			if nargin<4 || isempty(indz2)
				indz2 = 1:obj.NzT2;
			end
			if nargin<5 || isempty(indw)
				indw = obj.nw;
			end
			
			%	Get dimension sizes
			%	----------------------------------------------------------------
			if islogical(indz1)
				Nz1 = sum(indz1);
			else
				Nz1 = length(indz1);
			end
			
			%	Try loading quick method (only access subset of data in files -
			%	this has restrictions on valid indexing
			if exist('fast', 'var') && fast
				try
					I = cell2mat(reshape(cellfun(@(SD) ...
						SD.I(indth, indw, indz2), obj.ScanData(indz1), ...
						'UniformOutput', false), [1 1 1 Nz1]));
				catch
					fast = false;
				end
			end
			
			%	Load using 'slow' method - no restriction on how to index data
			%	subset, but loads all data from file
			if ~exist('fast', 'var') || ~fast
				I = cell2mat(reshape(cellfun(@(SD) subsref(SD.I, ...
					struct('type', '()', 'subs', {{indth, indw, indz2}})), ...
					obj.ScanData(indz1), 'UniformOutput', false), ...
					[1 1 1 Nz1]));
			end
			
			%	Permute data to I(th, z1, z2, w)
			I_r_z1_z2_w = permute(I, [1 4 3 2]);
		end

		%	====================================================================
		%		G E T   H A R M O N I C   L O C A T I O N S
		%	====================================================================
		function I = GetHarmonicLocations(obj, bg, Peaks, HarmWidth, I, Plot, ...
				varargin)
			
			%	Check peak parameters
			if ~exist('Peaks', 'var') || ~isstruct(Peaks)
				Peaks = struct;
			end
			if ~isfield(Peaks, 'Width')
				Peaks.Width = 50;
			end
			if ~isfield(Peaks, 'Height')
				Peaks.Height = 100;
			end
			if ~isfield(Peaks, 'Number')
				Peaks.Number = 8;
			end
			if ~isfield(Peaks, 'Indices')
				Peaks.Indices = 1:Peaks.Number;
			end
			
			%	Check harmonic width
			if ~exist('HarmWidth', 'var') || isempty(HarmWidth)
				HarmWidth = 5;
			end
			
			%	Load data if necessary
			%	----------------------------------------------------------------
			if ~exist('I', 'var') || isempty(I)
				I = zeros(obj.Nth, obj.Nw, 'uint32');
				
				%	Average measured data
				for nz=1:obj.NzT1
					I = I + sum(cast(obj.ScanData{nz}.I, 'like', I), ...
						3, 'native');
				end
			end
			
			%	Integrate over divergence
			if exist('bg', 'var') && ~isempty(bg)
				Iw = sum(obj.BG(I, bg, true), 1);
			else
				Iw = sum(I, 1);
			end

			%	Find central frequency of harmonics
			%		Should double check this for each data set!
			nw0 = sort(FindPeaks(Iw, Peaks.Width, Peaks.Height, ...
				Peaks.Number));
			Peaks.Indices(Peaks.Indices>length(nw0)) = [];
			obj.nw0 = nw0(Peaks.Indices);
			
			%	Find central divergence of harmonics
			[Imx, nth0] = max(mean(reshape( ...
				I(:, obj.nw0+(-HarmWidth:HarmWidth)), ...
				size(I, 1), obj.Nq, []), 3)); %#ok<*PROPLC>
			obj.nth0 = round(Weighted_Mean(nth0, Imx));
			
			%	Plot results
			%	----------------------------------------------------------------
			if exist('Plot', 'var') && Plot
				%	Setup figure
				fig = figure(gcf);
				clf(fig);
				fig.Position(3:4) = round(1e3*[1 1.2*size(I, 1)/size(I, 2)]);

				%	Plot harmonics
				imagesc(log10(double(obj.BG(I, ...
					FindArg(varargin, 'Background', 2.4e4)))));
				
				%	Annotate plot
				axis image;
				title(sprintf(['Harmonic positions: <\\theta_0> = %.3gpxl,' ...
					'\\sigma_{\\theta_0} = %.3gpxl'], obj.nth0, std(nth0)));
				xlabel('Frequency [pxl]');
				ylabel('Divergence [pxl]');
				colormap(gca, jet(2^8));

				%	Add positions
				line(obj.nw0, nth0, 'LineStyle', 'none', 'Marker', 'o', ...
					'Color', 'w', 'MarkerSize', 10, 'LineWidth', 2)
				line(obj.nw([1 end]), obj.nth0*[1 1], 'LineStyle', '--', ...
					'Color', 'w', 'LineWidth', 2);
			end
			
			
		end
		
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (Set / Get)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		
		%	====================================================================
		%		G E T   F O L D E R
		%	====================================================================
		function Folder = get.Folder(obj)
			%	Use base and raster folders only
			Folder = fullfile(obj.BaseFolder, obj.RasterFolder);
			if exist(Folder, 'dir')
				return;
			end
			
			%	Apply date between base and raster folders
			DateDir = regexpi(obj.RasterFolder, '\d{4}-\d{2}-\d{2}', ...
				'match', 'once');
			Folder = fullfile(obj.BaseFolder, DateDir, obj.RasterFolder);
			if ~exist(Folder, 'dir')
				error('Could not find folder');
			end
		end
		
		%	====================================================================
		%		G E T   I N F O
		%	====================================================================
		function Info = get.CamInfo(obj)
			str = GetFilenames(obj.Folder, '*CameraInfo.mat');
			Info = load(fullfile(obj.Folder, str{1}));
		end
		
		function Info = get.ScanInfo(obj)
			str = GetFilenames(obj.Folder, '*ScanInfo.mat');
			Info = load(fullfile(obj.Folder, str{1}));
		end

		
		%	====================================================================
		%		G E T   T A R G E T   P O S I T I O N S
		%	====================================================================
		function z = get.zT1(obj)
			z = Column(obj.ScanInfo.Positions{2}{1});
		end
		
		function z = get.zT2(obj)
			z = Column(obj.ScanInfo.Positions{1}{1});
		end
		
		function Nz = get.NzT1(obj)
			Nz =  obj.ScanInfo.Steps{2};
		end
		
		function Nz = get.NzT2(obj)
			Nz =  obj.ScanInfo.Steps{1};
		end
		
		%	====================================================================
		%		G E T   I M A G E   A X E S
		%	====================================================================
		function Nw = get.Nw(obj)
			Nw =  double(obj.CamInfo.Width);
		end
		
		function Nth = get.Nth(obj)
			Nth =  double(obj.CamInfo.Height);
		end
		
		function nw = get.nw(obj)
			nw = Column(1:obj.Nw);
		end
		
		function nth = get.nth(obj)
			nth = Column(1:obj.Nth);
		end
		
		%	====================================================================
		%		M I S C E L L A N E O U S
		%	====================================================================
		function Nq = get.Nq(obj)
			Nq = length(obj.nw0);
		end
		
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (Static)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Static)
		%	====================================================================
		%		B A C K G R O U N D   R E M O V A L
		%	====================================================================
		function I = BG(I, bg, floor)
			I = (I-bg);
			
			%	This is not necessary if input is unsigned integer (underflow is
			%	automatically truncated, unless Matlab's defauls settings have
			%	been modified).
			if nargin>2 && ~isempty(floor) && floor
				I(I<0) = cast(0, 'like', I);
			end
		end
		
		%	====================================================================
		%		C A P P I N G   /   B A C K G R O U N D   R E M O V A L
		%	====================================================================
		function I = Cap(I, cap)
			I = I.*cast(I<=cap, 'like', I) ...
				+ cast(cap, 'like', I)*cast(I>cap, 'like', I);
		end
	end
end