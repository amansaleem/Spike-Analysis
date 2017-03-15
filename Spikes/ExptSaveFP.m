function success = ExptSaveFP(exptfp,overwrite_flag)
% Saves the exptfp data in a mat file
%
% success = ExptSaveFP saves the exptfp in global PICK
%
% success = ExptSaveFP(exptfp) lets you specify the exptfp
% success = ExptSaveFP(exptfp,overwrite_flag) lets you overwrite an
% existing file (default is false);
%
% 2011-07 Neel Dhruv based on ExptSave.m by MC

global DIRS PICK;

if nargin < 1
    exptfp = [];
    if isfield(PICK,'exptfp');
        exptfp = PICK.exptfp;
    end
end
if nargin < 2
    overwrite_flag = false;
end

success = false; % start pessimistic
    
if isempty(exptfp)
    return
end

if ~isfield(exptfp,'Source')
    msgbox('No field "source"!!');
    return
end

AnimalDir = fullfile(DIRS.data,exptfp.animal);
if ~isdir(AnimalDir), mkdir(AnimalDir); end

SeriesDir = fullfile(AnimalDir,num2str(exptfp.iseries));
if ~isdir(SeriesDir), mkdir(SeriesDir); end

ExpDir = fullfile(SeriesDir,num2str(exptfp.iexp));
if ~isdir(ExpDir), mkdir(ExpDir); end

FileName = sprintf( 'Expt_%s.mat', exptfp.Source );

if exist(fullfile(ExpDir, FileName),'file') && ~overwrite_flag
    msgbox('File exists!!');
    return
end

% option v6 disables compression, thus saving time
fprintf('Saving exptfp file... ');
save( fullfile(ExpDir, FileName), 'exptfp', '-v6' );
fprintf('done\n');

success = true; 


