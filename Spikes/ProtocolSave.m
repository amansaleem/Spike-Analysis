function success = ProtocolSave( Protocol, quiet )
% ProtocolSave saves a Protocol in the appropriate directory
%
% success = ProtocolSave( Protocol )
%
% success = ProtocolSave( Protocol, 'quiet' ) suppresses any text
% output (important e.g. when called from a Visual Basic program).
%
% CAREFUL: it overwrites!
%
% 2006-08 Matteo Carandini

global DIRS


if nargin<2
    quiet = 'loud';
end

if ~isfield(DIRS,'data')
    error('No field "data" in DIRS');
end

success = 0;

animal = Protocol.animal;
iseries = Protocol.iseries;
iexp = Protocol.iexp;

if ~isdir(fullfile(DIRS.data,animal))
    mkdir(fullfile(DIRS.data,animal));
end

if ~isdir(fullfile(DIRS.data,animal,num2str(iseries)))
    mkdir(fullfile(DIRS.data,animal,num2str(iseries)));
end

if ~isdir(fullfile(DIRS.data,animal,num2str(iseries),num2str(iexp)))
    mkdir(fullfile(DIRS.data,animal,num2str(iseries),num2str(iexp)));
end

FileName = fullfile(DIRS.data,animal,num2str(iseries),num2str(iexp),'Protocol.mat');

% if exist( FileName, 'file' )
%     errordlg(sprintf('File %s exists!',FileName));
%     return
% end

if ~strcmp(quiet,'quiet')
    fprintf('--> Writing %s\n',FileName);
end

save(FileName,'Protocol');

success = 1;