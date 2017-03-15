function noise = ExptEstimateNoise(expt,ichan)
% ExptEstimateNoise
% 
% noise = ExptEstimateNoise(expt)
% noise = ExptEstimateNoise(expt,ichan)
% estimates the noise as described by Equation 7 of "Spike detection
% and sorting algorithms for polytrodes" by Blanche et al
%
% Example: 
% expt = ExptLoad('catz025',2,1);
% noise = ExptEstimateNoise(expt);
%
% 2009-02 Matteo Carandini

global PICK
global DIRS

if nargin < 1
    expt = PICK.expt;
end

if nargin < 2
    noise = zeros(expt.nchans,1);
    for ichan = 1:expt.nchans
       noise(ichan) = ExptEstimateNoise(expt,ichan);
    end
    return
end

if nargin == 2
    chans = ChanLoad(DIRS.spikes,expt.animal,expt.iseries,expt.iexp);
    chan = chans(ichan);
    expt = ExptFilter(expt,chan);
    
    WinDur = 0.01; % ms
    WinNSamples = floor(expt.samplerate*WinDur);
    mm = zeros(expt.nstim,expt.nrepeats);
    for is = 1:expt.nstim
        for ir = 1:expt.nrepeats
            TotNSamples = length(expt.data{ichan}{is,ir});
            nchunks = floor(TotNSamples/WinNSamples);
            ChunkMedians = zeros(nchunks,1);
            for ichunk = 1:nchunks
                ii = (1:WinNSamples)+(ichunk-1)*WinNSamples;
                ChunkMedians(ichunk) = median(abs(expt.data{ichan}{is,ir}(ii)));
            end
            mm(is,ir) = median(ChunkMedians);
        end
    end
    noise = median(mm(:))/0.6745;
    return
end


