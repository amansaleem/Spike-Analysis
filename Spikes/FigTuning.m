function fig = FigTuning()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.

load FigTuning

h0 = figure('Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'FileName','G:\Code\Spikes\FigTuning.m', ...
	'MenuBar','none', ...
	'NumberTitle','off', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperType','A4', ...
	'PaperUnits','points', ...
	'Position',[300 200 200 400], ...
	'Tag','Fig1', ...
	'ToolBar','none');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat1, ...
	'ListboxTop',0, ...
	'Position',[11.25 14.25 45 24.75], ...
	'Style','edit', ...
	'Tag','txtMaxResp');
if nargout > 0, fig = h0; end
