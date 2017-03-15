function [meanModel, bins, cellList, maxPos] = getPopNormVector(model, cellList, reorder)

idx = 1;
if nargin<3
    reorder = 0;
end
if nargin<2
    cellList = [];
    resort = 1;
else
    resort = 0;
end

maxPos = zeros(1,length(model.performance));
for icell = 1:length(model.performance)
    [~,maxPos(icell)] = max(nanmean(model.model.tuning(icell).respModel,1)*model.sampleRate);
end

if resort
    for icell = find(model.performance>0)
        if icell< 186
            meanModel(idx,:) = nanmean(model.model.tuning(icell).respModel,1)*model.sampleRate;
            cellList = [cellList icell];
            idx = idx + 1;
        end
    end
elseif reorder
    nlist = find(model.performance>0);
    badList = [];
    for icell = 1:length(cellList)
        if sum(cellList(icell)==nlist)==0
            badList=[badList icell];
        end
    end
    cellList(badList) = [];
    for icell = cellList
        if icell< 186
            meanModel(idx,:) = nanmean(model.model.tuning(icell).respModel,1)*model.sampleRate;
%             cellList = [cellList icell];
            idx = idx + 1;
        end
    end
else
    for icell = cellList
        if icell< 186
            meanModel(idx,:) = nanmean(model.model.tuning(icell).respModel,1)*model.sampleRate;
            %         cellList = [cellList icell];
            idx = idx + 1;
        end
    end
end

bins = model.bins;

if resort
    lowFiring = find(max(meanModel')<1);
    meanModel(lowFiring,:) = [];
    cellList(lowFiring) = [];
    
    
    [~,max_idx] = max(meanModel');
    [~,sort_idx] = sort(max_idx);
    
    meanModel = meanModel(sort_idx,:);
    cellList  = cellList(sort_idx);
end

meanModel = meanModel./(ones(50,1)*sum(meanModel'))';

% meanModel = (meanModel - (ones(50,1)*min(meanModel'))');
% meanModel = (meanModel ./ (ones(50,1)*max(meanModel'))');