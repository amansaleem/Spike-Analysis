function myscreen = ScreenLogLoad(animal,iseries,iexp,rateflag)
% ScreenLogLoad load the screenlog of an experiment
%
% myscreen = ScreenLogLoad( ExptTag ) loads the screenlog of the experiment
% specified by ExptTag.animal, ExptTag.iseries, ExptTag.iexp.
%
% myscreen = ScreenLogLoad(animal,iseries,iexp) does the same job.
%
% myscreen = ScreenLogLoad(animal,iseries,iexp,rateflag) lets you specify 
% rateflag: - 'real' (default) loads the real, physical framerate
%           - 'stim' loads the wrong framerate, used to generate stimuli
%
% 2003-02 VM made it
% 2006-05 MC first argument can be ExptTag with appropriate fields.
% 2007-04 MC changed so the imposition of 124 Hz applies only up to catz071
% 2007-05 MC because of bug in last edit 124 Hz imposition was never applied...
% 2009-02 AZ (UNDONE) motherdir is now two directories above DIRS.data in order to find the correct screendir - needed because \Cat and \Mouse designation
% 2010-07 ND added check for screen logs to be in DIRS.stiminfo - needed because of new stimbox saves screen logs there instead of in old location
% 2010-09 ND added a new hack to get a RealFrameRate field in for the shared rig setup
% 2011-03 ND hacked a way to allow the directory to be more general than
% the animal name, for example the stim info for M110301_BALL is in the
% folder M110301
% 2011-06 ND for some reason screen information now comes in as
% myScreenInfo rather than what it used to come in as (myscreen). Made it
% work with either new stuff or old stuff.
% 2012-04 ND made sure could load new screen info rather than just the dummy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    New screen information does not allow assigning a field to
%    RealFrameRate as used to be possible. Could make code not backwards
%    compatible. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% why is line 95 hard coded to be assuming 120 Hz? we don't use crt's anymore.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global DIRS;

if nargin == 1
    ExptTag = animal; % the first argument was ExptInfo
    animal  = ExptTag.animal;
    iseries = ExptTag.iseries;
    iexp 	= ExptTag.iexp;
end

if nargin < 4 || isempty(rateflag)
   rateflag = 'real';
end

gotit = 0;
if ((strcmpi(animal,'catz') && str2double(animal(end-2:end)) > 12) || ~strcmpi(animal,'catz'))
   motherdir = fileparts(DIRS.data);
   logname = sprintf('%s/%s_%i_%i.mat',upper(animal),lower(animal),iseries,iexp);
   screendir = fullfile(motherdir,'screen logs',upper(animal));
   if ~exist(screendir,'dir')
      warning(['Directory ' screendir ' does not exist']);
   else
      fprintf(1,'Loaded screen log file %s in directory %s\n',logname,screendir);
      %previousdir=pwd;
      %cd(screendir);
      screenfile = fullfile(motherdir,'screen logs',logname);
      gotit = exist(screenfile,'file');
      if gotit
         load(screenfile);
      end
      %cd(previousdir);
   end
end

if gotit~=2
    if ((strcmpi(animal,'catz') && str2double(animal(end-2:end)) > 12) || ~strcmpi(animal,'catz'))
        motherdir = DIRS.stimInfo;
        if isempty(strfind(animal,'_'))
            logname = sprintf('%s/%s_%i_%i.mat',upper(animal),lower(animal),iseries,iexp);
            screendir = fullfile(motherdir,upper(animal));
        else
            logname = sprintf('%s/%s_%i_%i.mat',upper(animal(1:strfind(animal,'_')-1)),lower(animal),iseries,iexp);
            screendir = fullfile(motherdir,upper(animal(1:strfind(animal,'_')-1)));
        end
%         screendir = fullfile(motherdir,upper(animal(1:strfind(animal,'_')-1)));
        if exist(screendir,'dir')~=7
            warning(['Directory ' screendir ' does not exist']);
        else
            fprintf(1,'Loaded screen log file %s in directory %s\n',logname,screendir);
            screenfile = fullfile(motherdir,[upper(logname(1:end-4)), '.mat']);
            gotit = exist(screenfile,'file');
            if gotit==2
                load(screenfile);
            end
        end
    end
end

if gotit ~= 2
   myscreen = [];
   wtext1 = sprintf('----->WARNING:Wasn`t able to load the screen log of %s %i %i',animal,iseries,iexp);
   wtext2 = sprintf('----->            will load a FICTICIOUS screen log         ');
   disp(wtext1);
   disp(wtext2);
end

%% Override the FrameRate

% it turns out that when going at about 125 Hz on the Sony Multiscan 500PS
% the refresh rate is always 124.8537 Hz, no matter what the screenlog says
% so let's impose this value.

% true when we were using the Mac (up to catz071):
if ( (strcmpi( animal(1:4),'catz') && str2double(animal(end-2:end)) < 71) )
    RealFrameRate = 124.8537;
else
    %RealFrameRate = myscreen.FrameRate;
    RealFrameRate = 120.3115;
end

% sometimes screen information in myscreen structure and sometimes in
% myScreenInfo ScreenInfo structure
if ~exist('myscreen','var')
    myscreen = myScreenInfo;
    % try catch bit commented out by ND 2012_04_12 - was always just
    % loading the dummy ScreenInfo and then moving on. 
%     try
%         myscreen = ScreenInfo;
%     catch
%         myscreen = myScreenInfo;
%     end
end

if ~exist('myscreen','var')
    % ND 2012_04_12 - I don't think this if statement is ever entered
    myscreen = myScreenInfo;
end

if ~isempty(myscreen) && strcmp(myscreen.MonitorType,'Multiscan 500PS')

   % The FrameRate used to generate the stimuli
   StimFrameRate = myscreen.FrameRate;

   switch rateflag
      case 'stim'
         myscreen.RealFrameRate = RealFrameRate;
         myscreen.StimFrameRate = StimFrameRate;

      case 'real'
         myscreen.FrameRate     = RealFrameRate;
         myscreen.RealFrameRate = RealFrameRate;
         myscreen.StimFrameRate = StimFrameRate;
   end
   
elseif ~isempty(myscreen) && strncmp(myscreen.MonitorType,'HannsG HW191D',length('HannsG HW191D')) % hack introduced by ND

   % The actual framerate of the monitor
   % problem assigning something called RealFrameRate in new type of screen
   % information. commenting out for now. could cause problems later.
   
%    myscreen.RealFrameRate = myscreen.FrameRate;
   
elseif isempty(myscreen)
    
   % default when we were using the Mac (up to catz071)
%    myscreen.WhichScreen = 1;
%    myscreen.PixelDepth = 8;
%    myscreen.windowPtr = NaN;
%    myscreen.ScreenRect = [0 0 640 480];
%    myscreen.Xmax = 640;
%    myscreen.Ymax = 480;
%    myscreen.FrameRate = RealFrameRate;
%    myscreen.MonitorType = 'Fictititious';
%    myscreen.MonitorSize = 39;
%    myscreen.CalibrationDir = 'None';
%    myscreen.PixelSize = 0.0609;
%    myscreen.Dist = 59;
%    myscreen.RealFrameRate = RealFrameRate;
%    myscreen.StimFrameRate = RealFrameRate;
   
   myscreen.WhichScreen = 1;
   myscreen.PixelDepth = 8;
   myscreen.windowPtr = NaN;
   myscreen.ScreenRect = [0 0 800 600];
   myscreen.Xmax = 800;
   myscreen.Ymax = 600;
   myscreen.FrameRate = 120.3115;
   myscreen.MonitorType = 'Fictitious';
   myscreen.MonitorSize = 39;
   myscreen.CalibrationDir = 'None';
   myscreen.PixelSize = 0.0488;
   myscreen.Dist = 10; % very small distance
   myscreen.RealFrameRate = 120.3115;
   myscreen.StimFrameRate = 120.3115;

end

% The old (beautiful) code:
% if strcmp( myscreen.MonitorType, 'Multiscan 500PS') & abs(myscreen.FrameRate-124)<2
%     myscreen.FrameRate = 124.8537;
% end


return

