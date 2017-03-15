function UnitSave(unit,spikedir)
% UnitSave saves a Unit data structure
%
% UnitSave(unit,spikedir)
%
% where spikedir is usually DIRS.spikes
% part of Spikes
%
% 2000-09 MC
% 2008-02 MC updated by making the directory (may be unnecessary)
% 2009-02 AZ Instead of overwriting an existing file, UnitSave now archives
% the old one in an 'archive' directory, with a timestamp: -yyyymmddHHMM
% 2010-06 MC added saving units with traces -- no archiving for those

if isfield(unit,'datatype') && ~isempty(unit.datatype)
    datatype = unit.datatype;
else
    if strcmp(unit.icell,'traces')
        datatype = 'traces';
    else
        datatype = 'spiketimes';
    end
end

filename = UnitGetFilename( unit.animal, unit.iseries, unit.iexp, unit.ichan, unit.icell );

TheDir = fullfile(spikedir,unit.animal,num2str(unit.iseries));
if ~isdir(TheDir), mkdir(TheDir); end


switch datatype
    case 'spiketimes'
        fprintf(1,'\tSaving spike data for channel %d, cell %d...',unit.ichan, unit.icell);
        % AZ 2009-02: if sorted unit already exists, archive old file with
        % timestamp, e.g., 'CATZ085\6\CATZ085_s06_e09_c1015_u005.mat'
        % becomes    'CATZ085\archive\CATZ085_s06_e09_c1015_u005-200902091328.mat'
        xst = exist([spikedir filesep filename], 'file');  % test if unit already exists
        if xst
            UnitArchive(unit,spikedir,filename);
        end
    case 'traces'
        fprintf(1,'\tSaving trace data for channel %d...',unit.ichan);
    otherwise
        error('Huh??');
end


% a = ver('MATLAB');
% MatlabVersion = str2num(a.Version)

a = version;
MatlabVersion = str2double(a(1));

if MatlabVersion < 7
    save(fullfile(spikedir,filename),'unit')
else
    save(fullfile(spikedir,filename),'unit','-v6')
end

fprintf(1,'done\n')
