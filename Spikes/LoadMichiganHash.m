function [hash, fs] = LoadMichiganHash(prot)
% Loads the high-pass filtered signals for Michigan probe recordings - no
% spike sorting just the hash
%
% [hash, fs] = LoadMichiganFPs(prot) returns a cell array of high-pass filtered
% signals and their sample rate in Hz (fs).
%
% hash is a cell array with dimensions [nstim,nrepeats]. Each cell has
% dimensions [nz,nsamples], where nsamples can vary with stimulus and repeat.
%
% 2007-09 Matteo Carandini
% 2009-11 ND made from LoadMichiganFPs

global DIRS

hash = cell( prot.nstim, prot.nrepeats );

for istim = 1:prot.nstim
    fprintf('Stim %d: ', istim);

    for irepeat = 1:prot.nrepeats
        fprintf('.');
        filename = sprintf('%s_%d_%d_%d_%d-Michigan.mat',prot.animal, prot.iseries, prot.iexp, irepeat, istim);
        datafile = fullfile(DIRS.data,prot.animal,int2str(prot.iseries),int2str(prot.iexp),filename);
        load(datafile,'stimresps');
        hash{istim,irepeat} = stimresps.data;
    end
        fprintf('\n');
end

fs = stimresps.samplerate;
