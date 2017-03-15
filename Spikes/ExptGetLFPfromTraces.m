function [expt,fps] = ExptGetLFPfromTraces(expt,desired_lfp_samplerate)
% 
% ExptsGetLFPfromTraces takes the raw data from expts.data and resamples to
% create lowpass filtered and downsampled copies of the data that isolates
% the LFP signal. Does this using the resample function in Matlab
%
% expts = ExptGetLFPfromTraces(expts)
% 
% expt.lfp.data is a cell array with dimensions [nstim,nrepeats]. Each cell has
%      dimensions [nchans,nsamples], where nsamples can vary with stimulus and repeat.
% 
% expt.lfp.samplerate is the new samplerate appropriate for the lfp's.
% 
% expt.lfp.filtercoeffs contains the filter coefficients that were used by
%      the resample function.
% 
% 2010-06 ND created
% 2010-09 ND let PICK be defined by expt data handed in
% 2011-07 ND changed samplerate of FPs from 500 Hz to 400 Hz
% 2011-07 ND actually allowed input argument to choose desired lfp
%            samplerate (default is 400 Hz).
% 2011-07 ND removed change from 2010-09 - don't really need PICK, just use
%            it in to define expt if expt is not already an input argument.

global DIRS
global PICK

% no input argument, so assume PICK defined and get data from there.
if nargin < 1
    if isempty(PICK.expt)
        error('<ExptGetLFPfromTraces> Could not find global variable PICK or no data in PICK.expt');
    end
    expt = PICK.expt;
end
% ND removed this
% % ND added because sometimes might call this function without using the
% % ExptPickFcn routines, in which case PICK is not defined.
% if isempty(PICK)
%     PICK.expt = expt;
% end

if nargin < 2
    desired_lfp_samplerate = 400;
end

%% get some information about nchans, nstims, nrepeats

nchans = size(expt.data,1);
[nstims,nrepeats] = size(expt.data{1});

% get original samplerate
orig_samplerate = expt.samplerate;

% dowsample from orig_samplerate to desired_lfp_samplerate
downsample_factor = floor(orig_samplerate/desired_lfp_samplerate);
true_lfp_samplerate = orig_samplerate/downsample_factor;

% resample for each channel and each stim and each repeat
fprintf('Computing LFP''s ');
expt.lfp.data = cell(nstims,nrepeats);
for istim = 1:nstims
    fprintf('.');
    for irepeat = 1:nrepeats
        for ichan = 1:nchans
            if isempty(expt.data{ichan}{istim,irepeat})
                continue;
            end
            [resample_output,filtercoeffs] = resample(double(expt.data{ichan}{istim,irepeat}),1,downsample_factor);
            expt.lfp.data{istim,irepeat}(ichan,:) = resample_output;
        end
    end
end
fprintf(' done.\n');

expt.lfp.samplerate = true_lfp_samplerate;
expt.lfp.filtercoeffs = filtercoeffs;
fps = expt.lfp;

return

%% if you want to see how it does, compare the following two plots

% this_stim = ceil(nstims.*rand(1,1));
% this_repeat = ceil(nrepeats.*rand(1,1));
% this_chan = ceil(nchans.*rand(1,1));
% 
% figure;
% ax(1) = subplot(221); plot(expt.data{this_chan}{this_stim,this_repeat},'b');
% [n_orig,bins] = hist(expt.data{this_chan}{this_stim,this_repeat});
% subplot(223); bar(bins,n_orig/sum(n_orig),'FaceColor','w','EdgeColor','b');
% ax(2) = subplot(222); plot(expt.lfp.data{this_stim,this_repeat}(this_chan,:),'r');
% [n_filtered] = hist(expt.lfp.data{this_stim,this_repeat}(this_chan,:),bins);
% subplot(224); bar(bins,n_filtered/sum(n_filtered),'FaceColor','w','EdgeColor','r');
% matchy(ax);
