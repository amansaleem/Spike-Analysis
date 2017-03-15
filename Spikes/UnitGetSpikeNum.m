function [SpikeNum, nrpts, nstim] = UnitGetSpikeNum( unit, analysisSpan)

data = unit.spiketimes;
nstim = size(data,1);
nrpts	= size(data,2);
SpikeNum=zeros(nstim,nrpts);

for istim=1:nstim
    for irpt=1:nrpts%nrpts
        SpikeNum(istim,irpt)=length(find(unit.spiketimes{istim,irpt}>=analysisSpan(1)...
            & unit.spiketimes{istim,irpt}<=analysisSpan(2)));
    end
end
