function [PVO, binsA, binsB] = calculatePopVectorOverlap(modelA, modelB)

idx = 1;
for icell = find(modelB.performance>0 & modelA.performance>0)
    if icell< 186
        meanModelA(idx,:) = nanmean(modelA.model.tuning(icell).respModel,1)*modelA.sampleRate;
        meanModelB(idx,:) = nanmean(modelB.model.tuning(icell).respModel,1)*modelB.sampleRate;
        idx = idx + 1;
    end
end

binsA = modelA.bins;
binsB = modelB.bins;

lowFiring = find(max(meanModelA')<1 | max(meanModelB')<1);
meanModelA(lowFiring,:) = [];
meanModelB(lowFiring,:) = [];

[~, max_idx] = max(meanModelA');
[~,sort_idx] = sort(max_idx);

meanModelA = meanModelA(sort_idx,:);
meanModelB = meanModelB(sort_idx,:);

meanModelA = meanModelA./(ones(50,1)*sum(meanModelA'))';
meanModelB = meanModelB./(ones(50,1)*sum(meanModelB'))';

for ibinA = 1:size(meanModelA,2)
    for ibinB = 1:size(meanModelB,2)
        
        % From Rossavaad,... Mehta, Science, 2013
        PVO(ibinA, ibinB) = (meanModelA(:,ibinA)'*meanModelB(:,ibinB)) ./ ...
            (meanModelA(:,ibinA)'*meanModelA(:,ibinA))*(meanModelB(:,ibinB)'*meanModelB(:,ibinB));
        
    end
end