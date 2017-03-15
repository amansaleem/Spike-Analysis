function MoveDataToServer(animal)
% MoveDataToServer moved data to the server
%
% MoveDataToServer(animal) moves data from experiment machines to the
% server in the directories defined by SetDefaultDirs
%
% Needs to be improved! For now moves only Screen Logs and data by
% Multispike
%
% 2009-09 Matteo Carandini and Tatsuo Sato
% 2009-10 input the names of the host PCs. by TS
% 2009-11 ND added backup of Michigan tanks
% 2010-06 ND added backup of Cerebus data
% 2010-06 ND added ability to enter 'none' for a host and then it wouldn't
% back that one up

if nargin < 1
    error('Need an animal name');
end


%% the name of the host PCs  (2009-10 by TS)
NameStimHost=[];
while isempty(NameStimHost)
    NameStimHost=input('stim host PC  >>','s');
end
NameDaqHost=[];
while isempty(NameDaqHost)
    NameDaqHost=input('data acq host PC  >>','s');
end
NameTankHost=[];
while isempty(NameTankHost)
    NameTankHost=input('michigan tank host PC  >>','s');
end
NameCerebusHost=[];
while isempty(NameCerebusHost)
    NameCerebusHost=input('cerebus host PC  >>','s');
end

SetDefaultDirs

%% the screen logs

if ~strcmp(NameStimHost,'none')
    % SourceDir = fullfile('\\zap\Screen Logs',animal);
    SourceDir = fullfile(strcat('\\',NameStimHost,'\Screen Logs'),animal);% 2009-10 by TS
    DestinDir = fullfile('\\Zserver\data\screen logs',animal);



    if ~isdir(DestinDir)
        status = mkdir(DestinDir);
    end

    if ~isdir(SourceDir) || ~isdir(DestinDir)
        MsgBox('Problem!!!');
    end

    DosString = sprintf('xcopy "%s\\*.*" "%s" /e ', SourceDir, DestinDir );
    [status,result] = dos( DosString, '-echo' );
end

%% the raw data files from Multispike

if ~strcmp(NameDaqHost,'none')
    % SourceDir = fullfile('\\zazzera\data\data',animal);
    SourceDir = fullfile(strcat('\\',NameDaqHost,'\data\data'),animal);% 2009-10 by TS
    DestinDir = fullfile(DIRS.data,animal);

    if ~isdir(DestinDir)
        status = mkdir(DestinDir);
    end

    if ~isdir(SourceDir) || ~isdir(DestinDir)
        MsgBox('Problem!!!');
    end

    DosString = sprintf('xcopy "%s\\*.*" "%s" /e ', SourceDir, DestinDir );
    [status,result] = dos( DosString, '-echo' );
end

%% the Michigan data

if ~strcmp(NameTankHost,'none')
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
end

%% the Blackrock data

if ~strcmp(NameCerebusHost,'none')
    % SourceDir = fullfile('\\zanna\cerebus_data',animal);
    SourceDir = fullfile(strcat('\\',NameCerebusHost,'\cerebus_data'),animal);% 2010-06 by ND
    DestinDir = fullfile(DIRS.Cerebus,animal);

    if ~isdir(DestinDir)
        status = mkdir(DestinDir);
    end

    if ~isdir(SourceDir) || ~isdir(DestinDir)
        MsgBox('Problem!!!');
    end

    DosString = sprintf('xcopy "%s\\*.*" "%s" /e ', SourceDir, DestinDir );
    [status,result] = dos( DosString, '-echo' );
end
