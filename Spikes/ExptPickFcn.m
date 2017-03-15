function varargout = ExptPickFcn(animal,iseries,iexp)
%
% EXPTPICKFCN (formerly SPIKES) defines a few environment variables and launches the Graphical User Interfaces
%

% 2009-01 SK 
% Small change to allow the user to store sorted spikes in
% a user-defined directory: SetDefaultDirs is only called if DIRS does not
% exist.
% 2009-02 AZ
% Turned spikes.m into a function.
% 2009-03 AZ
% Renamed spikes and split into function and script versions. With this FUNCTION 
% version, the user can specify a specific animal, and optionally, series,
% and experiment to open with:
% >> ExptPickFcn(animal)
% >> ExptPickFcn(animal,iseries)
% >> ExptPickFcn(animal,iseries,iexp)
% Can either declare no output variables as above, or half, or all of them:
% >> [DEMO,DIRS,PICK] = ExptPickFcn
% >> [DEMO,DIRS,PICK,animal,iexp,iseries] = ExptPickFcn(animal,iseries,iexp)

% addpath('c:/toolbox\matteobox');
% addpath('c:/users/matteo/matlab/toolbox5\matteobox');
%
% 1) Go to control panel > system > advanced > environment variables 
% 2) Choose New, and make two new variables, one called DATADIR, and one called 
% SPIKEDIR, with values appropriate to your configuration (in my case, Z:\trodes and S:\). No quote marks.
% 3) Restart Matlab if it was running. 

% to make sure that analringach works:
global DEMO
DEMO = 1;

global DIRS PICK

% -------- default data and spikes directories
% % AZ 2009-02: SetDefaultDirs if any DIRS field is nonexistent (UNDONE)
if ~exist('DIRS', 'var') || isempty(DIRS) %|| ~isfield(DIRS,'data') || ...
%         ~isfield(DIRS,'spikes')   || ~isfield(DIRS,'camera')   || ...
%         ~isfield(DIRS,'xfiles')   || ~isfield(DIRS,'michigan') || ...
%         ~isfield(DIRS,'Cerebus')  || ~isfield(DIRS,'stiminfo')
    SetDefaultDirs;
end

% % -------- override with environment variables, if defined ----------------
% 
% [problemflag, dirname]=system('echo %DATADIR%');
% dirname = dirname(find(dirname+1-1~=10)); % remove character '10' inserted by Windows
% if ~strcmp(dirname,'%DATADIR%'), DIRS.data = dirname; end
% 
% [problemflag, dirname]=system('echo %SPIKEDIR%');
% dirname = dirname(find(dirname+1-1~=10)); % remove character '10' inserted by Windows
% if ~strcmp(dirname,'%SPIKEDIR%'), DIRS.spikes = dirname; end
% 
% [problemflag, dirname]=system('echo %CAMERADIR%');
% dirname = dirname(find(dirname+1-1~=10)); % remove character '10' inserted by Windows
% if ~strcmp(dirname,'%CAMERADIR%'), DIRS.camera = dirname; end

%--------------------------------------------------------------------------

while ~exist(DIRS.data,'dir')
   answer = inputdlg('Cannot open DATA directory. Please enter a new one:',...
      'Spikes',1,{DIRS.data});
   if isempty(answer), return; end
   DIRS.data = answer{1};
end

while ~exist(DIRS.spikes,'dir')
   answer = inputdlg('Cannot open SPIKES directory. Please enter a new one:',...
      'Spikes',1,{DIRS.spikes});
   if isempty(answer), return; end
   DIRS.spikes = answer{1};
end

while ~exist(DIRS.camera,'dir')
   answer = inputdlg('Cannot open CAMERA directory. Please enter a new one:',...
      'Spikes',1,{DIRS.spikes});
   if isempty(answer), return; end
   DIRS.camera = answer{1};
end

PICK.expt = [];
PICK.chans = [];
% AZ 2010-05-08 Needed to bypass pbAnimal
if isfield(PICK,'exptinfos')
   PICK = rmfield(PICK,'exptinfos');
end

% AZ 2009-02-25: User-defined picks
switch nargin
   case 0
      % do nothing
   case 1
      PICK.animal  = animal;
   case 2
      PICK.animal  = animal;
      PICK.iseries = iseries;
   case 3
      PICK.animal  = animal;
      PICK.iseries = iseries;
      PICK.iexp    = iexp;
   otherwise
      warning('Too many input arguments!');
end

%------- Launch the GUIs

% PICK.picker = FigPicker;

PICK.picker = openfig('FigPick.fig','reuse');
set(PICK.picker,'HandleVisibility','callback');

PICK.spikesorter = FigSpikeSorter;
set(PICK.spikesorter,'HandleVisibility','callback');
set(PICK.spikesorter,'Visible','Off','WindowStyle','Normal','RendererMode','auto','resize','off'); % 'CloseRequestFcn',@DenyInteractiveCloseRequest);



% AZ 2009-02-25 If user specified animal, pretend we clicked the 'Animal' button
if isfield(PICK,'animal') && ~isempty(PICK.animal)
   FigPicker_callbacks pbAnimal;
end

% AZ 2009-03-06 don't return anything if no output arguments specified
switch nargout
   case 0
      % do nothing
   case 3
      varargout = {DEMO,DIRS,PICK};
   case 6
      varargout = {DEMO,DIRS,PICK,animal,iexp,iseries};
   otherwise
      % do nothing
      disp([num2str(nargout),' output variables passed.  Try again.'])
end

end % end function ExptPickFcn(...)