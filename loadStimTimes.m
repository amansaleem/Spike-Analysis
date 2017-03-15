function [screenTimes, stimOn, stimOff] = loadStimTimes(animal, iseries, iexp)

SetDefaultDirs
global pepNEV
nevDir = fullfile(DIRS.Cerebus,animal,num2str(iseries));

KlustersDir = fullfile(DIRS.multichanspikes,animal,num2str(iseries));

% nevSamplingRateInKHZ = 30;
period = 1;%nevSamplingRateInKHZ*1000 / samplingRate;

nevFileName = [nevDir filesep num2str(iexp) filesep animal '_' num2str(iseries) '_' num2str(iexp) '.nev'];
if exist(nevFileName, 'file')
    nevopen(nevFileName);
    screenTimes  = pepNEV.sync.timestamps;
    nevclose;
else
    disp('Trying to load by repeat...')
    screenTimes  = [];
    irep = 1;
    nevFileName = [nevDir filesep num2str(iexp) filesep animal ...
        '_' num2str(iseries) '_' num2str(iexp) '_' num2str(irep) '.nev'];
    while exist(nevFileName, 'file')
        nevopen(nevFileName);
        screenTimes  = [screenTimes  pepNEV.sync.timestamps];
        nevclose;
        irep = irep + 1;
        nevFileName = [nevDir filesep num2str(iexp) filesep animal ...
            '_' num2str(iseries) '_' num2str(iexp) '_' num2str(irep) '.nev'];
    end
end
% disp('Removed condition to remove stimuli within 100 ms');
screenTimes(find(diff(screenTimes)<100) + 1) = [];

stimOn = screenTimes(1:2:end);
stimOff = screenTimes(2:2:end);
