function createKwikUnits(animal, iseries, iexp, igroup, addInfo)

SetDefaultDirs
% DIRS.spikes = [DIRS.spikes filesep 'Klustered'];

chans = getKwikSpikes(animal, iseries, iexp, igroup, addInfo);
% [stimTimes,p] = createStimMatrix(animal, iseries, iexp);
p = ProtocolLoad(animal, iseries, iexp);
[~, stimOnTimes, stimOffTimes] = loadStimTimes(animal, iseries, iexp);

numClus = length(chans);
uu(numClus) = UnitCreate;
period = chans(1).sampleRate;

for clusID = 1:numClus
    uu(clusID).animal       = animal;
    uu(clusID).iseries      = iseries;
    uu(clusID).iexp         = iexp;
    uu(clusID).id           = chans(clusID).id;
    uu(clusID).ichan        = igroup;
    uu(clusID).icell        = chans(clusID).icell; %clusID
    uu(clusID).datatype     = 'spiketimes';
    uu(clusID).nstim        = p.nstim;
    uu(clusID).nrepeats     = p.nrepeats;
    uu(clusID).prestimdurs  = [];
    uu(clusID).waveformMean = [];
    uu(clusID).waveformStd  = [];
    uu(clusID).isolDist     = [];
    uu(clusID).lRatios      = [];
    uu(clusID).prototype    = [];
    uu(clusID).sampledur    = [];
    uu(clusID).source       = 'kwik file';
    for istim = 1 : p.nstim
        for irep =  1 : p.nrepeats
            t0 = stimOnTimes(p.seqnums(istim,irep));
            t1 = stimOffTimes(p.seqnums(istim,irep));
            if p.seqnums(istim,irep) > 1
                tmin1 = stimOffTimes(p.seqnums(istim,irep)-1);
            else
                tmin1 = 0;
            end
            
            % convert to units of samples
            t0 = (t0/period);
            t1 = (t1/period);
            tmin1 = (tmin1/period);
            
            % find spike times belonging to this trial
          
            currCellSpikes= chans(clusID).spiketimes;
            spikeTimes    = currCellSpikes(currCellSpikes>=t0 & currCellSpikes<=t1 );
            preSpikeTimes = currCellSpikes(currCellSpikes>=tmin1 & currCellSpikes<=t0);
            
            spikeTimes    = spikeTimes - t0; % time = 0 is now the start of the stimulus
            preSpikeTimes = preSpikeTimes - t0;
            
            uu(clusID).spiketimes{istim,irep} = spikeTimes';
            uu(clusID).prespiketimes{istim,irep} = preSpikeTimes';
            uu(clusID).stimdurs(istim,irep)     = t1-t0;
        end
        clear spikeTimes preSpikeTimes waveForms
    end
end
for clusID = 1:numClus
    UnitSave_Imperial(uu(clusID),[DIRS.spikes filesep 'Klustered' ]);
end