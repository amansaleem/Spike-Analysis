function mixModels = calculateMixtureModels(es, t)

theta = 0:pi/8:pi-0.1;

for iangle = 1:length(theta)
    stim = sin(theta(iangle))*es.smthBallSpd(t) + cos(theta(iangle))*es.smthTrajSpd(t);
    mixModels.angle(iangle) = theta(iangle);
    
    mixModels.data(iangle).map = oneDimMap;
    mixModels.data(iangle).map = mixModels.data(iangle).map.trainSpikeMap(stim, es.spikeTrain(t,:), 150);
    
    mixModels.performance(iangle,:) = mixModels.data(iangle).map.performance;
    display(['Completed angle ' num2str(iangle) ' of ' num2str(length(theta))]);
end
