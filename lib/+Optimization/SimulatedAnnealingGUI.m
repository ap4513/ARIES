classdef SimulatedAnnealingGUI < uiobjects.HardwareGUI
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	properties
	end

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	methods

		%	====================================================================
		
		function obj = SimulatedAnnealingGUI(varargin)
			obj@uiobjects.HardwareGUI(varargin);
% 			obj.CreateGUI;
		end
		
		%	====================================================================
		
	end %	methods
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (Access=protected)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	methods (Access=protected)

		%	====================================================================

		function CreateGUI(obj, varargin)
			%	Call superclass method (generates main GUI alignment controls)
			CreateGUI@uiobjects.HardwareGUI(obj);

			%	Load XML file
			%	----------------------------------------------------------------
			xmlfile = obj.FindArg(varargin, 'XMLFile');
			if isempty(xmlfile)
				xmlfile = fullfile(...
					obj.GetFileLocation('\+Optimization'), ...
					'SimulatedAnnealingGUI.xml');
			end
			
			%	Add control panels to interface
			%	----------------------------------------------------------------
			obj.AddChildren(obj.LoadGUI(xmlfile), ...
				obj.Handles.ControlPanelAlignmentVBox);

		end


	end %	methods (Access=protected)
end %	classdef
