function [model1 model2 model3 spred1 spred2 spred3] = get2Dmap(zua, varA, varB, nGrid, binsA, binsB, CVO, iter, sampFreq, win, win2)

% [model] = get1Dmap(zua, variable, nGrid, bins, CVO, iter, sampFreq, win)

X_train(:,1) = varA(CVO.train{iter});
X_test(:,1)  = varA(CVO.test{iter});
X_cv(:,1)    = varA(CVO.cv{iter});

X_train(:,2) = varB(CVO.train{iter});
X_test(:,2)  = varB(CVO.test{iter});
X_cv(:,2)    = varB(CVO.cv{iter});

pos_on_map_cv      = nGrid(1)*(X_cv(:,2)-1) + X_cv(:,1);
pos_on_map      = nGrid(1)*(X_test(:,2)-1) + X_test(:,1);

if nargout<2
    processAll = 0;
else
    processAll = 1;
end
% for icell = 1:numCells
if nargin<9
    sampFreq = 200;
    win = 10;
    smth = 0;
else
    smth = 1;
end

if nargin<11
    win2 = win;
end

if smth
    zua_train = smthInTime(zua(CVO.train{iter}), sampFreq, win);
    zua_test  = smthInTime(zua(CVO.test{iter}), sampFreq, win);
    zua_cv    = smthInTime(zua(CVO.cv{iter}), sampFreq, win);
else
    zua_train = zua(CVO.train{iter});
    zua_test  = zua(CVO.test{iter});
    zua_cv    = zua(CVO.cv{iter});
end
test = zua(CVO.test{iter});


%get spike count map
scMap = full(sparse(X_train(:,1), X_train(:,2), (zua_train), nGrid(1), nGrid(2)));
% get occupancy map
occMap = full(sparse(X_train(:,1), X_train(:,2), 1, nGrid(1), nGrid(2)));

if processAll
    % These are the marginals for the additive model
    scMapA = full(sparse(X_train(:,1), 1, (zua_train), nGrid(1), 1));
    occMapA = full(sparse(X_train(:,1), 1, 1, nGrid(1), 1));
    scMapB = full(sparse(X_train(:,2), 1, (zua_train), nGrid(2), 1));
    occMapB = full(sparse(X_train(:,2), 1, 1, nGrid(2), 1));
end

n1 = 10;
grids1 = 1:n1;
n2 = 10;
grids2 = 1:n2;

EV1 = zeros(length(grids1),length(grids2));
if processAll
    EV2 = zeros(length(grids1),length(grids2));
    EV3 = zeros(length(grids1),length(grids2));
end
for win1=1:length(grids1)
    for win2=1:length(grids2)
        % Full joint model
        FRMap = special_smooth_2d(scMap, [1./grids1(win1) 1./grids2(win2)], binsA, binsB, [n1 n2])...
            ./special_smooth_2d(occMap, [1./grids1(win1) 1./grids2(win2)], binsA, binsB,[n1 n2]);
        fullFRvector    = reshape(FRMap, [],1);
        pred1 = fullFRvector(pos_on_map_cv);
        if smth
            spred1  = smthInTime(pred1, sampFreq, win);
        else
            spred1 = pred1;
        end
        EV1(win1,win2) = calCrossValExpVar(zua_train, zua_cv, spred1);
        
        if processAll
            % Multiplying SVD
            trash = nanmean(FRMap(:));
            if isnan(trash)
                trash =0;
            end
            FRMap(isnan(FRMap)) = trash;
            [~,~,~,FRMap_svd] = MakeSeparable(FRMap,0);
            fullFR_svd_vector    = reshape(FRMap_svd, [],1);
            pred2 = fullFR_svd_vector(pos_on_map_cv);
            if smth
                spred2  = smthInTime(pred2, sampFreq, win);
            else
                spred2 = pred2;
            end
            EV2(win1,win2) = calCrossValExpVar(zua_train, zua_cv, spred2);
            
            % Multiplying marginals
            FRMap1 = special_smooth_1d(scMapA, 1./grids1(win1), binsA, n1)...
                ./special_smooth_1d(occMapA, 1./grids1(win1), binsA, n1);
            FRMap2 = special_smooth_1d(scMapB, 1./grids2(win2), binsB, n2)...
                ./special_smooth_1d(occMapB, 1./grids2(win2), binsB, n2);
            pred3 = FRMap1(X_cv(:,1)).*FRMap2(X_cv(:,2));
            if smth
                spred3  = smthInTime(pred3, sampFreq, win);
            else
                spred3 = pred3;
            end
            EV3(win1,win2) = calCrossValExpVar(zua_train, zua_cv, spred3);
        end
    end
end

EV1 = round(100*EV1);
if processAll
    EV2 = round(100*EV2);
    EV3 = round(100*EV3);
end

[idxX1, idxY1]= find(EV1 == max(max(EV1)),1,'first');
if processAll
    [idxX2, idxY2]= find(EV2 == max(max(EV2)),1,'first');
    [idxX3, idxY3]= find(EV3 == max(max(EV3)),1,'first');
end

% Full joint model
if ~(isempty(idxX1) | isempty(idxY1) | isnan(idxX1) | isnan(idxY1))
    
    [model1.swin(1)] = grids1((idxX1));
    [model1.swin(2)] = grids2((idxY1));
    model1.bins{1} = binsA;
    model1.bins{2} = binsB;
    
    model1.tuning = special_smooth_2d(scMap, [1./model1.swin(1) 1./model1.swin(2)], binsA, binsB, [n1 n2])...
        ./special_smooth_2d(occMap, [1./model1.swin(1) 1./model1.swin(2)], binsA, binsB, [n1 n2]);
    %             pred  = model.tuning(X_test);
    fullFRvector    = reshape(model1.tuning, [],1);
    %                     pred  = fastsmooth(fullFRvector(pos_on_map),FRsmthwin,3,1);
    pred = fullFRvector(pos_on_map);
    if smth
        spred1  = smthInTime(pred, sampFreq, win);
    else
        spred1 = pred1;
    end
    [model1.EV model1.corr model1.L model1.Q model1.train_mean] = calCrossValExpVar(zua_train, zua_test, spred1, test, pred);
    
    %             [icell model.swin model.EV model.corr]
else
    model1.swin = nan;
    model1.tuning = nan;
    model1.EV = nan;
    model1.corr = nan;
    model1.L = nan;
    model1.Q = nan;
    model1.bins{1} = binsA;
    model1.bins{2} = binsB;
    
end

if processAll
    % Multiplying svd model
    if ~(isempty(idxX2) | isempty(idxY2) | isnan(idxX2) | isnan(idxY2))
        
        [model2.swin(1)] = grids1((idxX2));
        [model2.swin(2)] = grids2((idxY2));
        model2.bins{1} = binsA;
        model2.bins{2} = binsB;
        
        FRMap = special_smooth_2d(scMap, [1./model2.swin(1) 1./model2.swin(2)], binsA, binsB, [n1 n2])...
            ./special_smooth_2d(occMap, [1./model2.swin(1) 1./model2.swin(2)], binsA, binsB, [n1 n2]);
        
        trash = nanmean(FRMap(:));
        if isnan(trash)
            trash =0;
        end
        FRMap(isnan(FRMap)) = trash;
        [~,~,~,model2.tuning] = MakeSeparable(FRMap,0);
        
        fullFRvector    = reshape(model2.tuning, [],1);
        %                     pred  = fastsmooth(fullFRvector(pos_on_map),FRsmthwin,3,1);
        pred = fullFRvector(pos_on_map);
        if smth
            spred2  = smthInTime(pred, sampFreq, win);
        else
            spred2 = pred2;
        end
        [model2.EV model2.corr model2.L model2.Q] = calCrossValExpVar(zua_train, zua_test, spred2, test, pred);
    else
        
        model2.swin = nan;
        model2.tuning = nan;
        model2.EV = nan;
        model2.corr = nan;
        model2.L = nan;
        model2.Q = nan;
        model2.bins{1} = binsA;
        model2.bins{2} = binsB;
        
    end
    
    % Multiplying marginals
    if ~(isempty(idxX3) | isempty(idxY3) | isnan(idxX3) | isnan(idxY3))
        
        [model3.swin(1)] = grids1((idxX3));
        [model3.swin(2)] = grids2((idxY3));
        model3.bins{1} = binsA;
        model3.bins{2} = binsB;
        
        model3.tuning{1} = special_smooth_1d(scMapA, 1./model3.swin(1), binsA, n1)...
            ./special_smooth_1d(occMapA, 1./model3.swin(1), binsA, n1);
        model3.tuning{2} = special_smooth_1d(scMapB, 1./model3.swin(2), binsB, n2)...
            ./special_smooth_1d(occMapB, 1./model3.swin(2), binsB, n2);
        
        pred = model3.tuning{1}(X_test(:,1)).*model3.tuning{2}(X_test(:,2));
        if smth
            spred3  = smthInTime(pred, sampFreq, win);
        else
            spred3 = pred3;
        end
        [model3.EV model3.corr model3.L model3.Q] = calCrossValExpVar(zua_train, zua_test, spred3, test, pred);
    else
        model3.swin = nan;
        model3.tuning = nan;
        model3.EV = nan;
        model3.corr = nan;
        model3.L = nan;
        model3.Q = nan;
        model3.bins{1} = binsA;
        model3.bins{2} = binsB;
    end
end
