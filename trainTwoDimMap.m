function [obj, Prediction, V] = trainTwoDimMap(obj, V, Y, smth_win);

if isempty(obj.kfold)
    obj.kfold = 5;
end

if isempty(obj.train_mean)
    calTrainMean = 1;
else
    calTrainMean = 0;
end

if sum(isnan(V))>0 | sum(isnan(Y(:)))>0
    display('WARNING!!! Nans in the data: making a temp fix');
    t = ones(size(V,1));
    t(isnan(sum(V,2))) = 0;
    t(isnan(sum(Y,2))) = 0;
    V = V(t>0,:);
    Y = Y(t>0,:);
end

if sum(V<0)>0
    display('WARNING!!! -ve variable data present: ignorning them');
    t = ones(size(V,1));
    t(V(:,1)<0) = 0;
    t(V(:,2)<0) = 0;
    V = V(t>0,:);
    Y = Y(t>0,:);
end
if isempty(obj.CVO)
    obj.CVO = crossValPartition(ones(1,length(V)),obj.kfold);
end

if obj.kfold == 1
    obj.CVO.kfold = 1;
    obj.CVO.train{1}= ones(1,length(V));
    obj.CVO.cv{1}   = obj.CVO.train{1};
    obj.CVO.test{1} = obj.CVO.train{1};
end
for icell = 1:size(Y,2)
    Y(:,icell) = smthInTime(Y(:,icell), obj.sampleRate, smth_win);
end

obj.performance = zeros(1,obj.CVO.kfold);
Prediction = [];

% Bayes decoder needs discretization
[V(:,1) obj.bins]   = normalise1var(V(:,1), obj.numBins);
[V(:,2) obj.binsB]  = normalise1var(V(:,2), obj.numBinsB);

for iter = 1:obj.CVO.kfold
    
    Ytrain  = Y(obj.CVO.train{iter},:);
    Ycv     = Y(obj.CVO.cv{iter},:);
    Ytest   = Y(obj.CVO.test{iter},:);
    
    Vtrain  = V(obj.CVO.train{iter},:);
    Vcv     = V(obj.CVO.cv{iter},:);
    Vtest   = V(obj.CVO.test{iter},:);
    
    if calTrainMean
        train_mean = mean(Vtrain,1);
    else
        train_mean = obj.train_mean;
    end
    
    % Getting the 1D map for each neuron (calculating the place fields)
    % ...the main training component
          display(['Processing iter: ' num2str(iter)]);
    for icell = 1:size(Y,2)
%         display(['Processing cell: ' num2str(icell)]);
        %[model] = get1Dmap(Y(:,icell), V', obj.numBins, obj.bins, obj.CVO, iter, obj.sampleRate, obj.smth_win);
        [model] = get2Dmap(Y(:,icell), V(:,1), V(:,2), [obj.numBins obj.numBinsB], obj.bins, obj.binsB, obj.CVO, iter, obj.sampleRate, obj.smth_win, obj.smth_win);
        
        if isnan(model.swin)
            obj.model.tuning(icell).respModel(iter,:,:) = NaN*ones(obj.numBins, obj.numBinsB);
        else    
            obj.model.tuning(icell).respModel(iter,:,:) = model.tuning;
        end
        obj.model.EV(iter,icell) = model.EV;
    end
    obj.model.EV(iter,obj.model.EV(iter,:)<0) = 0;
    
    %%
    Prediction = [];
end

obj.performance = mean(obj.model.EV,1);