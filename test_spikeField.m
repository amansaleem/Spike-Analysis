[es] = getVRspikes('M130920_BALL',1025,103);
t2 = ~isnan(es.ballspeed);
spd = zeros(size(es.ballspeed));
spd(t2) = smthInTime(es.ballspeed(t2), 60, 500);
t = es.traj>0 & spd>3 & es.contrast>0 & es.roomLength == 1;
k = bayesDecoder;
k.numBins = 100;
[k, prediction, X, Posterior] = trainDecoder(k, es.traj(t), es.spikeTrain(t,:));
for icell = (find( mean(k.model.EV)>0.05))
    subplot(211)
    plot(es.traj(t),es.sampleTimes(t), '.k', 'MarkerSize',2);
    hold on;
    plot(es.traj(es.spikeTrain(:,icell) & t),es.sampleTimes(es.spikeTrain(:,icell) & t), '.r', 'MarkerSize',15);
    hold off
    subplot(212)
    for iter = 1:length(k.model.trained)
        plot(k.model.trained(iter).respModel_orig(icell,:)*60);
        hold on;
    end
    hold off
    title(['Cell: ' num2str(icell) ', CellID: ' es.spikeIDs{icell} ' Mean Rate: ' num2str(nanmean(es.spikeTrain(:,icell))*60)])
    pause
end