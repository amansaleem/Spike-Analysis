function [fps, fs] = LoadMichiganFPs(prot,stim_ix)
% Loads the field potentials for Michigan probe recordings
%
% [fps, fs] = LoadMichiganFPs(prot) returns a cell array of field
% potentials (fps) and their sample rate in Hz (fs).
% 
% [fps, fs] = LoadMichiganFPs(prot,stim_ix) allows you to select just a
% particular set of the stims.
%
% fps is a cell array with dimensions [nstim,nrepeats]. Each cell has
% dimensions [nz,nsamples], where nsamples can vary with stimulus and repeat.
%
% 2007-09 Matteo Carandini
% 2010-06 ND added try statement so wouldn't break if a particular
%         stim/repeat were never converted. old version archived in etc folder
%         also added ability to select just some of the stims.

if nargin < 2, stim_ix = 1:prot.nstim; end

global DIRS

fps = cell( length(stim_ix), prot.nrepeats );

nstims = length(stim_ix);
for istim = 1:nstims
    this_stim = stim_ix(istim);
    fprintf('Stim %d: ', this_stim);

    for irepeat = 1:prot.nrepeats
        fprintf('.');
        filename = sprintf('%s_%d_%d_%d_%d-Michigan.mat',prot.animal, prot.iseries, prot.iexp, irepeat, this_stim);
        datafile = fullfile(DIRS.data,prot.animal,int2str(prot.iseries),int2str(prot.iexp),filename);
        try load(datafile,'lfp'); % added by ND
            fps{istim,irepeat} = lfp.data;
        catch
            warning('empty stim repeat. couldn''t load.');
            fps{istim,irepeat} = [];
        end
    end
    fprintf('\n');
end

fs = lfp.samplerate;
