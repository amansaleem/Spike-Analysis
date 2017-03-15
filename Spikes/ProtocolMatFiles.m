function [MatFileNames, AbsentMatFiles] = ProtocolMatFiles(protocol)
% ProtocolMatFiles returns a matrix of file names and gives a report
%
% ProtocolMatFiles(protocol) displays a report about absent mat files
%
% MatFileNames = ProtocolAbsentMatFiles(protocol) returns a cell array 
%
% [MatFileNames, AbsentMatFiles] = ProtocolMatFiles(protocol) returns also
% array of logicals
%
% 2009-05 Matteo Carandini took code from ProtocolLoad
% 2010-06 MC converted from ProtocolAbsentMatFiles, added file name output

global DIRS

fprintf('Testing for completeness of data files...');

MatFileNames = cell(protocol.nstim,protocol.nrepeats);
AbsentMatFiles = false(protocol.nstim,protocol.nrepeats);
for irepeat = 1:protocol.nrepeats
    for istim = 1:protocol.nstim
        filename = sprintf('%s_%d_%d_%d_%d.mat',protocol.animal, protocol.iseries, protocol.iexp, irepeat, istim);
        datafile = fullfile(DIRS.data,protocol.animal,int2str(protocol.iseries),int2str(protocol.iexp),filename);
        if exist(datafile,'file')
            MatFileNames{istim,irepeat} = datafile;
        else
            AbsentMatFiles(istim,irepeat) = true;
            fprintf(' (no %s)', filename);
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
