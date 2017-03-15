function MoveCerebusToServer(animal,NameCerebusHost,verbose_flag,series,experiment)
% MoveCerebusToServer moves blackrock data to the server
%
% MoveCerebusToServer(animal) moves data from the machine that collects the
% data via the blackrock hardware and moves to the server in the
% directories defined by SetDefaultDirs NameCerebusHost is an optional
% argument for the computer hosting the original raw cerebus data. default
% is empty --> asks for user input. verbose_flag is an optional argument.
% default == 1 --> be verbose. series and experiment are optional arguments
% for selecting particular runs to be copied. useful so that not always
% copying entire directory contents (default --> copy all contents).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% should we make the xcopy occur with the option \Y --> no confirmation prompt?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% 2010-06 ND made this from MoveMichigantoServer
% 2010-07 ND made optional arguments NameCerebusHost and verbose_flag and
%            optional argument for series and experiment

if nargin < 1
    error('Need an animal name');
end
if nargin < 2
    NameCerebusHost=[];
end
if nargin < 3 || isempty(verbose_flag)
    verbose_flag = 1;
end
if nargin < 4
    series = [];
end
if nargin < 5
    experiment = [];
end

while isempty(NameCerebusHost)
    NameCerebusHost=input('cerebus host PC  >>','s');
end

SetDefaultDirs

%% Move the Cerebus data

% why has zpopulation been hardcoded in here. this is not right.
% SourceDir = fullfile('\\zpopulation\cerebus_data',animal);

% this is old way of saving one file per experiment and just one animal directory
SourceDir = fullfile(strcat('\\',NameCerebusHost,'\cerebus_data'),animal);
DestinDir = fullfile(DIRS.Cerebus,animal);
% this is new way of saving per repeat and each series and experiment is a folder
% make sure series and experiment are strings
% if isnumeric(series), series = num2str(series); end
% if isnumeric(experiment), experiment = num2str(experiment); end
% SourceDir = fullfile(strcat('\\',NameCerebusHost,'\cerebus_data'),animal,series,experiment);
% DestinDir = fullfile(DIRS.Cerebus,animal,series,experiment);

if ~isdir(DestinDir)
    status = mkdir(DestinDir);
end

if ~isdir(SourceDir) || ~isdir(DestinDir)
    MsgBox('Problem!!!');
end

% this is old way of saving one file per experiment and just one animal directory
% make sure series and experiment are scalars
if ischar(series), series = str2num(series); end
if ischar(experiment), experiment = str2num(experiment); end
if ~isempty(series) && ~isempty(experiment)
    SeriesExptString = ['u' sprintf('%03d',series) '_' sprintf('%03d',experiment)];
elseif ~isempty(series)
    SeriesExptString = ['u' sprintf('%03d*',series)];
else
    SeriesExptString = ['*'];
end
% this is new way of saving per repeat and each series and experiment is a folder
% SeriesExptString = ['*'];

if verbose_flag
    DosString = sprintf('xcopy "%s\\%s.*" "%s" /e ', SourceDir, SeriesExptString, DestinDir );
else
    DosString = sprintf('xcopy "%s\\%s.*" "%s" /e /y', SourceDir, SeriesExptString, DestinDir );
end
[status,result] = dos( DosString, '-echo' );
