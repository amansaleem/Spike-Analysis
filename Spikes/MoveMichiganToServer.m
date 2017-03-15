function MoveMichiganToServer(animal)
% MoveMichiganToServer moves michigan data tanks to the server
%
% MoveMichiganToServer(animal) moves data from the machine that collects the raw michigan data
% and moves to the server in the directories defined by SetDefaultDirs
%
% 2009-11 ND made this from MoveDatatoServer

if nargin < 1
    error('Need an animal name');
end


NameTankHost=[];
while isempty(NameTankHost)
    NameTankHost=input('michigan tank host PC  >>','s');
end

SetDefaultDirs

%% Move the Michigan data

% SourceDir = fullfile('\\zlack\data\tanks',animal);
SourceDir = fullfile(strcat('\\',NameTankHost,'\data\tanks'),animal);% 2009-11 by ND
DestinDir = fullfile(DIRS.michigan,animal);

if ~isdir(DestinDir)
    status = mkdir(DestinDir);
end

if ~isdir(SourceDir) || ~isdir(DestinDir)
    MsgBox('Problem!!!');
end

DosString = sprintf('xcopy "%s\\*.*" "%s" /e ', SourceDir, DestinDir );
[status,result] = dos( DosString, '-echo' );
