function [models] = getModels(animal, iseries,...
    only_correct, smthWin, box_filt, shank_list, loadA, area)

%% Parameters that can be changed
test_normal = 1; % Learning/training on normal condition
plot_all = 0;
numSampleBins = 250;

display(['Smoothing window = ' num2str(smthWin)]);
if nargin<3
    only_correct = false;
end
if nargin<4
    smthWin = 10;
end
if nargin<5
    box_filt = 1;
end
if nargin<6
    shank_list = false;
end
if nargin<7
    loadA = false;
end
if nargin<8
    area = 'CA1';
end

if length(num2str(animal))<=5
    animal = ['M130' num2str(animal) '_BALL'];
end

switch animal
    case 'M130918_BALL'
        iseries = 1030;
        switch area
            case 'V1'
                cell_list = 190:321;
            case 'CA1'
                cell_list = 1:189;
        end
    case 'M130920_BALL'
        iseries = 1025;
        switch area
            case 'V1'
                cell_list = 153:249;
            case 'CA1'
                cell_list = 1:152;
        end
    case 'M130703_BALL'
        iseries = 809;
        cell_list = [];
    otherwise
        if nargin<2
            iseries = inputdlg('Enter the Series:','0');
            iseries = str2num(iseries{1});
        end
        cell_list = [];
end
delayT = 0;

%% Load the spiking and theta data
es = [];

if loadA%iseries==601 | iseries==604 | iseries==1025 | iseries==530 | iseries==531
    load(['Data' filesep 'es_' animal '_' num2str(iseries) '_thetaBinsA_new2']);
else
    load(['Data' filesep 'es_' animal '_' num2str(iseries) '_thetaBinsB_new2']);
end

% load(['Data' filesep 'es_' animal '_' num2str(iseries) '_norm']);

switch iseries
    case 1025
        display('Warning!!! Removing interneurons on 1025');
        interneurons = zeros(1,size(es.spikeTrain,2));
        interneurons([17, 18,  9,11,12, 44, 45, 41, 78, 102, 104, 129, 113, 140, 145]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    
    case 1030
        display('Warning!!! Removing interneurons on 1030');
        interneurons = zeros(1,size(es.spikeTrain,2));
        interneurons([15, 20, 30, 39, 71, 132]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    
    case 602
        display('Warning!!! Removing interneurons on 602');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([5 8 10 18 24 29 35 38 40 52 54 59 60 61 77 90 92 93 107 109 127]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    case 604
        display('Warning!!! Removing interneurons on 604');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([8 22 25 26 19 31 40 52 53 66 91 97 109]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    case 531
        display('Warning!!! Removing interneurons on 531');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([11 15 16 44 53 86 94 98]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    case 603
        display('Warning!!! Removing interneurons on 603');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([3 8 9 15 57 63 64 66 74 116 131 138 171 177 182 185 188]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    case 601
        display('Warning!!! Removing interneurons on 601');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([4 9 10 11 12 20 24 29 47 50 52 60 64 71 84]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    case 530
        display('Warning!!! Removing interneurons on 530');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([4 8 36 39 56 62 72 78 84 95 98 101 99 102 103 105 131 155 163 167 207 217 222 229 231 266 270]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
end

% spdDec = bayesDecoder;
% spdDec.CVO = [];
% temp1 = es.smthBallSpd;
% % temp1(temp1<5) = nan;
% [spdDec,a] = spdDec.trainDecoder(temp1, es.spikeTrain,0);
% a = [a' nan*ones(1,length(es.smthBallSpd)-length(a))]';
% es.spikeTrain = a;
% 
% es.orig.spikeTrain = es.orig.smthBallSpd;

display(['Processing: ' animal '_' num2str(iseries)]);
%% Getting conditions
Posterior_all.animal = animal;
Posterior_all.iseries = iseries;

        base = es.traj>1 & es.contrast~=0 & ~isnan(es.traj)...
            & round(100*es.gain)/100==1 ...
            & round(100*es.roomLength)/100==1 ...
            & es.smthBallSpd>5 & es.trajspeed>=0;
        t      = es.traj~=0 & round(es.contrast*100)/100==0.6 & base;
        t_low  = es.traj~=0 & round(es.contrast*100)/100==0.18 & base;
        t_high = es.traj~=0 & round(es.contrast*100)/100==0.72 & base;
        t_gray = es.traj~=0 & es.contrast==0 & ~isnan(es.traj) & es.smthBallSpd>5;
        if strcmp(animal,'M130703_BALL');
            t      = es.contrast==0.75 & base;
            t_low  = es.contrast<0.75 & base;
            t_high = es.contrast>0.75 & base;
        end


spkRate = zeros(size(es.spikeTrain));
es.spikeTrain(es.traj==0,:) = nan; % to see if there is an effect of the gray screen.
for icell = 1:size(es.spikeTrain,2);
   spkRate(:,icell) = es.spikeTrain(:,icell);
end

%% Getting the models
if only_correct
    [model] = fitmodelTemp(es.traj(t & es.outcome==2), spkRate(t & es.outcome==2,:), smthWin);
    [model_low] = fitmodelTemp(es.traj(t_low & es.outcome==2), spkRate(t_low & es.outcome==2,:), smthWin);
    [model_high] = fitmodelTemp(es.traj(t_high & es.outcome==2), spkRate(t_high & es.outcome==2,:), smthWin);
else
    [model] = fitmodelTemp(es.traj(t), spkRate(t,:), smthWin);
    [model_low] = fitmodelTemp(es.traj(t_low), spkRate(t_low,:), smthWin);
    [model_high] = fitmodelTemp(es.traj(t_high), spkRate(t_high,:), smthWin);
end

models.norm = model;
models.low = model_low;
models.high = model_high;