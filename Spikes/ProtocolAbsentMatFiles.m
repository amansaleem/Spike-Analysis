
function AbsentMatFiles = ProtocolAbsentMatFiles(protocol)
% ProtocolAbsentMatFiles see which files are present and gives a report
%
% ProtocolAbsentMatFiles(protocol) displays the report
%
% AbsentMatFiles = ProtocolAbsentMatFiles(protocol) returns a matrix
%
% 2009-05 Matteo Carandini took code from ProtocolLoad
% 2010-06 MC added matrix output

global DIRS

fprintf('Testing for completeness of data files...');

AbsentMatFiles = false(protocol.nstim,protocol.nrepeats);
for irepeat = 1:protocol.nrepeats
    for istim = 1:protocol.nstim
        filename = sprintf('%s_%d_%d_%d_%d.mat',protocol.animal, protocol.iseries, protocol.iexp, irepeat, istim);
        datafile = fullfile(DIRS.data,protocol.animal,int2str(protocol.iseries),int2str(protocol.iexp),filename);
        if ~exist(datafile,'file')
            AbsentMatFiles(istim,irepeat) = true;
        end
    end
end

if any(AbsentMatFiles(:))
    for irepeat = 1:protocol.nrepeats
        nAbsentInThisRepeat = nnz(find(AbsentMatFiles(:,irepeat)));
        if nAbsentInThisRepeat>0
            fprintf('\n------> % d files missing in repeat %d\n',nAbsentInThisRepeat,irepeat);
        end
    end
else
    fprintf('done\n');
end
