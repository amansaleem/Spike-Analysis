function success = ExptSave(expt,lfp)
% Saves the expt data in a mat file
%
% success = ExptSave saves the expt in global PICK
%
% success = ExptSave(expt) lets you specify the expt
%
% 2010-06 Matteo Carandini
% 2011-07 MC upped saving option to v7.3 (slower but needed for large files)
% 2011-07 MC added the second argument, lfp (so it saves a second file)

global DIRS PICK;

if nargin < 1
    expt = [];
    if isfield(PICK,'expt');
        expt = PICK.expt;
    end
end

success = false; % start pessimistic
    
if isempty(expt)
    return
end

if ~isfield(expt,'Source')
    msgbox('No field "source"!!');
    return
end

AnimalDir = fullfile(DIRS.data,expt.animal);
if ~isdir(AnimalDir), mkdir(AnimalDir); end

SeriesDir = fullfile(AnimalDir,num2str(expt.iseries));
if ~isdir(SeriesDir), mkdir(SeriesDir); end

ExpDir = fullfile(SeriesDir,num2str(expt.iexp));
if ~isdir(ExpDir), mkdir(ExpDir); end

FileName = sprintf( 'Expt_%s.mat', expt.Source );
FileNameLFP = sprintf( 'Expt_%s_LFP.mat', expt.Source );

% if exist(fullfile(ExpDir, FileName),'file')
%     msgbox('File exists!! Not going to be saved');
%     return
% end

% MC 2011-12-21 added this information
foo = whos('expt');
MegaBytes = (foo.bytes)/2^20;
fprintf('The expt structure is %4.1f MBytes\n',MegaBytes);

% MC 2011-07-12 replaced option v6 (disables compression) 
% with option v7.3 (allows saving of larger files) 
fprintf('Saving expt file %s... ',fullfile(ExpDir, FileName));
save( fullfile(ExpDir, FileName), 'expt', '-v7.3' );
if nargin == 2 
    save( fullfile(ExpDir, FileNameLFP), 'lfp', '-v7.3' );
end

fprintf('done\n');

success = true; 


