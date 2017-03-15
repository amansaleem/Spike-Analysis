mGrid = 30;
pGrid = 30;

[mRate mBins] = normalise1var(modelInput(icell).meanrate, mGrid);
[pRate pBins] = normalise1var(modelInput(icell).poprate, pGrid);

zua = modelInput(icell).activity;
[sRate sBins] = normalise1var(zua, aGrid);

CVO = crossValPartition(ones(1,length(modelInput(icell).activity)),10);

for iter = 1:CVO.kfold
    [model(iexp).M(icell,iter) pred.M(icell,iter).pred] = get1Dmap(zua, mRate, mGrid, mBins, CVO, iter);
    
    [model(iexp).PM(icell,iter) model(iexp).PxM(icell,iter) model(iexp).PmM(icell,iter) ...
        pred.PM(icell,iter).pred pred.PxM(icell,iter).pred pred.PmM(icell,iter).pred] =...
        get2Dmap(zua, pRate, mRate, [pGrid mGrid], pBins, mBins, CVO, iter);
end