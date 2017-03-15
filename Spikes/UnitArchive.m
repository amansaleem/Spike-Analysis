function UnitArchive(unit,spikedir,filename)
% If sorted unit already exists, archive old file with
% timestamp, e.g., 'CATZ085\6\CATZ085_s06_e09_c1015_u005.mat'
% becomes    'CATZ085\archive\CATZ085_s06_e09_c1015_u005-200902091328.mat'
% Created by AZ 2009-02
% 2009-07    AZ made into function
% 2010-03-10 AZ Added series tag to discrim_pars files
%               Changed so timestamp is now(), rather than unit creation
%                 time. NOTE: time zones not supported

if nargin < 3
    filename = UnitGetFilename( unit.animal, unit.iseries, unit.iexp, unit.ichan, unit.icell );
end

old = load(fullfile(spikedir,filename)); % load existing unit
idx = find(filename==filesep);
old.filename = [unit.animal filesep 'archive' filesep filename(idx(end)+1:end)];
idx = find(old.filename=='.');

% AZ 2010-03-10: add -s01 if from series 1 (only if file is not a unit)
if strfind(filename,'discrim_pars')
   serIX = strfind(filename,filesep);
   iseries = str2double(filename(serIX(end-1)+1:serIX(end)-1));
   old.filename = sprintf('%s-s%02u%s',old.filename(1:idx(end)-1),iseries,...
      old.filename(idx(end):end));
end

% Timestamp: -yyyymmddHHMM
% AZ 2010-03-10: was: datestr(old.unit.timestamp,'yyyymmddHHMM')
old.filename = [old.filename(1:idx(end)-1) '-' ...
    datestr(now,'yyyymmddHHMM') old.filename(idx(end):end)];
if ~isdir([spikedir filesep unit.animal filesep 'archive'])
    mkdir([spikedir filesep unit.animal filesep 'archive']);
end
movefile(fullfile(spikedir,filename),fullfile(spikedir,old.filename));
fprintf('      Previously sorted unit\n      ''%s''\n      moved to\n      ''%s''\n',...
    filename,old.filename);