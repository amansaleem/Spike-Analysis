function [outputMatrix, p, cellIDs, outputTimes] = getStimSpiketimes(animal, iseries, iexp, shank_id, addInfo)

if nargin<4
    shank_id = 8;
end
if nargout>3
   getTimes = true;
else
   getTimes = false;
end
if nargin<5
    addInfo = '1';
end

% Load spike times and stim times
chans = getKwikSpikes(animal, iseries, iexp, shank_id, addInfo);
[stimTimes, p] = createStimMatrix(animal, iseries, iexp);

stimTimes.on = stimTimes.on./30000;
stimTimes.off = stimTimes.off./30000;

numStim = size(stimTimes.on,1);
numReps = size(stimTimes.on,2);
cellIDs = zeros(1,length(chans));

% stimTimes.on is istim x irep
for icell = 1:length(chans)
    st = chans(icell).spiketimes;
    cellIDs(icell) = chans(icell).icell;
    if length(st)<50
        continue
    end
    outputMatrix{icell} = zeros(numStim, numReps);
    for istim = 1:numStim
        for irep = 1:numReps
            startTime = stimTimes.on(istim,irep);
            endTime   = stimTimes.off(istim,irep);
            if getTimes
                outputTimes(icell, istim, irep).times = st(st>(startTime-0.5) & st<=(endTime+0.5)) - startTime;
            end
            outputMatrix{icell}(istim, irep) = sum(st>startTime & st<=endTime);
        end
    end
end
    
end