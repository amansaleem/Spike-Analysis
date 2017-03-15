function [stimTimes,p] = createStimMatrix(animal, iseries, iexp)

[~, stimOn, stimOff] = loadStimTimes(animal, iseries, iexp);

p = ProtocolLoad(animal, iseries, iexp);
[ss ii]=sort(p.seqnums(:));
seq=repmat([1:p.nstim]', p.nrepeats, 1);
seq=seq(ii);

stimTimes.on = zeros(p.nstim, p.nrepeats);
stimTimes.off= zeros(p.nstim, p.nrepeats);
    
for istim= 1:p.nstim
    stimTimes.on(istim,:)  = stimOn(seq==istim);
    stimTimes.off(istim,:) = stimOff(seq==istim);
end
