function [Posterior_all] = plot_some_decoding_new(animal, iseries, exp_list, type, only_correct, smthWin, box_filt, shank_list, area)

%% Parameters that can be changed
test_normal = 1; % Learning/training on normal condition
plot_all = 1;

display(['Smoothing window = ' num2str(smthWin)]);
if nargin<4
    type = 'contrast';
end
if nargin<5
    only_correct = false;
end
if nargin<6
    smthWin = 250;
end
if nargin<7
    box_filt = 0;
end
if nargin<8
    shank_list = false;
end
if nargin<9
    area = 'CA1';
end

if length(num2str(animal))<=5
    animal = ['M130' num2str(animal) '_BALL'];
end

if strcmp(animal,'M130918_BALL') && (nargin<2 || isempty(exp_list))
    exp_list = 103:105;
elseif strcmp(animal,'M130920_BALL') && (nargin<2 || isempty(exp_list))
    exp_list = 102:103;
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
% delayT = inputdlg('Enter delay frames:','s');
% delayT = str2num(delayT{1});
delayT = 0;

SetDirs;
es = VRLoadMultipleExpts(animal, iseries, exp_list,'SPIKES',[],shank_list);
if ~isempty(cell_list)
    es.spikeTrain = es.spikeTrain(:,cell_list);
end
es.spikeTrain = circshift(es.spikeTrain,[-delayT 0]);

%% Getting conditions

Posterior_all.animal = animal;
Posterior_all.iseries = iseries;
Posterior_all.exp_list = exp_list;

switch type
    case 'contrast'
        base = es.traj~=0 & es.contrast~=0 & ~isnan(es.traj) & es.gain==1 & es.roomLength==1 & es.smthBallSpd>5 & es.trajspeed>=0;
        t      = es.contrast==0.6 & base;
        t_low  = es.contrast==0.18 & base;
        t_high = es.contrast==0.72 & base;
        t_gray = es.traj~=0 & es.contrast==0 & ~isnan(es.traj) & es.smthBallSpd>5;
        if strcmp(animal,'M130703_BALL');
            t      = es.contrast==0.75 & base;
            t_low  = es.contrast<0.75 & base;
            t_high = es.contrast>0.75 & base;
        end
    case 'gain'
        base = es.traj~=0 & ~isnan(es.traj) & es.contrast==0.6 & es.roomLength==1 & es.smthBallSpd>5 & es.trajspeed>=0;
        t      = es.gain==1 & base;
        t_low  = es.gain<1 & base;
        t_high = es.gain>1 & base;
    case 'roomlength'
        base = es.traj~=0 & ~isnan(es.traj) & es.contrast==0.6 & es.gain==1 & es.smthBallSpd>5 & es.trajspeed>=0;
        t      = es.roomLength==1 & base;
        t_low  = es.roomLength<1 & base;
        t_high = es.roomLength>1 & base;
end

spkRate = zeros(size(es.spikeTrain));
for icell = 1:size(es.spikeTrain,2);
    if box_filt
        spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), 60, smthWin,'box');
    else
        spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), 60, smthWin); %*(smthWin*60./1000);
    end
end


%% Decoding all the conditions

dec = bayesDecoder;
dec.numBins = 50;

% if ~test_normal
dec_l = dec;
dec_h = dec;
% end
if only_correct
    dec2 = dec;
    [dec, ~, ~, ~] = ...
        dec.trainDecoder( es.traj(t & es.outcome==2), spkRate(t & es.outcome==2,:), 0);
    [~,  Posterior_norm] = dec.predictBayesDecoder(spkRate(t,:), 0,'best');
    [dec2, ~, X_norm, ~] = ...
        dec2.trainDecoder( es.traj(t), spkRate(t,:), 0);
else
    [dec, ~, X_norm, Posterior_norm] = ...
        dec.trainDecoder( es.traj(t), spkRate(t,:), 0);
end

if ~test_normal
    [dec_l, ~, X_low, Posterior_low] = ...
        dec_l.trainDecoder( es.traj(t_low), spkRate(t_low,:), 0);
    [dec_h, ~, X_high, Posterior_high] = ...
        dec_h.trainDecoder( es.traj(t_high), spkRate(t_high,:), 0);
    Posterior_high_orig =[];
    Posterior_low_orig  =[];
else
    [dec_l, ~, X_low_orig, Posterior_low_orig] = ...
        dec_l.trainDecoder( es.traj(t_low), spkRate(t_low,:), 0);
    [dec_h, ~, X_high_orig, Posterior_high_orig] = ...
        dec_h.trainDecoder( es.traj(t_high), spkRate(t_high,:), 0);
    [~,  Posterior_low] = dec.predictBayesDecoder(spkRate(t_low,:), 0,'best');
    [~,  Posterior_gray] = dec.predictBayesDecoder(spkRate(t_gray,:), 0,'best');
    [~,  Posterior_high] = dec.predictBayesDecoder(spkRate(t_high,:), 0,'best');
end

if test_normal
    [X_low, bins_low] = normalise1var(es.traj(t_low), dec.numBins);
    [X_high,bins_high] = normalise1var(es.traj(t_high), dec.numBins);
    X_high_orig = X_high(1:size(Posterior_high_orig,1));
    X_low_orig  = X_low(1:size(Posterior_low_orig,1));
else
    X_high = X_high(1:size(Posterior_high,1));
    bins_low = dec_l.bins;
    X_low  = X_low(1:size(Posterior_low,1));
    bins_high = dec_h.bins;
end
X_norm = X_norm(1:size(Posterior_norm,1));

timePs = find(t);
timePs = timePs(1:size(Posterior_norm,1));
t = false(size(t));
t(timePs) = true;

timePs = find(t_low);
timePs = timePs(1:size(Posterior_low,1));
t_low = false(size(t_low));
t_low(timePs) = true;

timePs = find(t_high);
timePs = timePs(1:size(Posterior_high,1));
t_high = false(size(t_high));
t_high(timePs) = true;

%% Defining the outputs
Posterior_all.Posterior_norm = Posterior_norm;
Posterior_all.Posterior_high = Posterior_high;
Posterior_all.Posterior_low  = Posterior_low ;

Posterior_all.Posterior_gray  = Posterior_gray;
Posterior_all.only_correct = only_correct;

Posterior_all.X_norm = X_norm;
Posterior_all.X_high = X_high;
Posterior_all.X_low  = X_low;

Posterior_all.data      = es;
Posterior_all.decoder   = dec;
Posterior_all.decoder_low   = dec_l;
Posterior_all.decoder_high   = dec_h;
Posterior_all.t_norm    = t;
Posterior_all.t_low     = t_low;
Posterior_all.t_high    = t_high;

Posterior_all.meanrates.low = mean(es.spikeTrain(t_low,:),1);
Posterior_all.meanrates.norm = mean(es.spikeTrain(t,:),1);
Posterior_all.meanrates.high = mean(es.spikeTrain(t_high,:),1);

Posterior_all.Posterior_high_orig = Posterior_high_orig;
Posterior_all.Posterior_low_orig  = Posterior_low_orig;
Posterior_all.X_high_orig = X_high_orig;
Posterior_all.X_low_orig  = X_low_orig;

%% Finding the trials as correct, early, miss and late
trialEnds = [find(diff(es.trialID)>=1)-15];
temp = zeros(size(es.traj));
temp(trialEnds) = 1;
trialEnds = temp;

outcome.complete        = (es.trajPercent(trialEnds>0)>80) | es.outcome(trialEnds>0)~=1;
outcome.correctTrials   = es.trialID(trialEnds & es.outcome==2);
outcome.earlyTrials     = es.trialID(trialEnds & (es.outcome==0 & es.trajPercent<80));
outcome.misslateTrials  = es.trialID(trialEnds  & es.outcome==0 & es.trajPercent>80);
outcome.lowContrast     = es.trialID(trialEnds & t_low);
outcome.normContrast    = es.trialID(trialEnds & t);
outcome.highContrast    = es.trialID(trialEnds & t_high);

es.outcome(es.outcome==1) = NaN;

% Splitting miss and late trials
outcome.missTrials = [];
outcome.lateTrials = [];
for tIdx = outcome.misslateTrials'
    % Does the animal lick after the reward position?
    if sum(es.lick((es.trialID==tIdx) & (es.trajPercent>70)))>0
        outcome.lateTrials = [outcome.lateTrials tIdx];
    else
        outcome.missTrials = [outcome.missTrials tIdx];
    end
end

% Early
for tIdx = [outcome.earlyTrials']
    es.outcome(es.trialID==tIdx) = 1;
end

% Correct
for tIdx = [outcome.correctTrials']
    es.outcome(es.trialID==tIdx) = 2;
end
% Late
for tIdx = [outcome.lateTrials]
    es.outcome(es.trialID==tIdx) = 3;
end
% Miss
for tIdx = [outcome.missTrials]
    es.outcome(es.trialID==tIdx) = 4;
end

%% Getting the confusion matrixes
for n = 1:dec.numBins
    %
    timePts = (X_norm == n);
    meanPost_norm(:,n) = nanmean(Posterior_norm(timePts, :),1);
    %
    timePts = (X_high == n);
    meanPost_high(:,n) = nanmean(Posterior_high(timePts, :),1);
    %
    timePts = (X_low == n);
    meanPost_low(:,n) = nanmean(Posterior_low(timePts, :),1);
    %
    timePts = (X_high_orig == n);
    meanPost_high_orig(:,n) = nanmean(Posterior_high_orig(timePts, :),1);
    %
    timePts = (X_low_orig == n);
    meanPost_low_orig(:,n) = nanmean(Posterior_low_orig(timePts, :),1);
    %
    timePts = (X_norm == n & (es.outcome(t) ~= 2));
    meanPost_norm_pass(:,n) = nanmean(Posterior_norm(timePts, :),1);
    %
    timePts = (X_norm == n & (es.outcome(t) == 2));
    meanPost_norm_actv(:,n) = nanmean(Posterior_norm(timePts, :),1);
    %
    timePts = (X_low == n & (es.outcome(t_low) == 2));
    meanPost_low_actv(:,n) = nanmean(Posterior_low(timePts, :),1);
    %
    timePts = (X_low == n & (es.outcome(t_low) ~= 2));
    meanPost_low_pass(:,n) = nanmean(Posterior_low(timePts, :),1);
    %
    timePts = (X_high == n & (es.outcome(t_high) ~= 2));
    meanPost_high_pass(:,n) = nanmean(Posterior_high(timePts, :),1);
    %
    timePts = (X_high == n & (es.outcome(t_high) == 2));
    meanPost_high_actv(:,n) = nanmean(Posterior_high(timePts, :),1);
end

outcome.allContrast = [outcome.lowContrast' outcome.normContrast' outcome.highContrast']';

%% Getting the Early, Correct, Late and Miss trials
for n = 1:dec.numBins
    % Early
    % % Low
    timePts1 = (X_low == n & (es.outcome(t_low) == 1));
    meanPost.low.early(:,n) = nanmean(Posterior_low(timePts1, :),1);
    % % Norm
    timePts2 = (X_norm == n & (es.outcome(t) == 1));
    meanPost.norm.early(:,n) = nanmean(Posterior_norm(timePts2, :),1);
    % % High
    timePts3 = (X_high == n & (es.outcome(t_high) == 1));
    meanPost.high.early(:,n) = nanmean(Posterior_high(timePts3, :),1);
    % % All
    meanPost.all.early(:,n) = nanmean([Posterior_low(timePts1, :)' Posterior_norm(timePts2, :)' Posterior_high(timePts3, :)']',1);
    
    % Correct
    % % Low
    timePts1 = (X_low == n & (es.outcome(t_low) == 2));
    meanPost.low.correct(:,n) = nanmean(Posterior_low(timePts1, :),1);
    % % Norm
    timePts2 = (X_norm == n & (es.outcome(t) == 2));
    meanPost.norm.correct(:,n) = nanmean(Posterior_norm(timePts2, :),1);
    % % High
    timePts3 = (X_high == n & (es.outcome(t_high) == 2));
    meanPost.high.correct(:,n) = nanmean(Posterior_high(timePts3, :),1);
    % % All
    meanPost.all.correct(:,n) = nanmean([Posterior_low(timePts1, :)' Posterior_norm(timePts2, :)' Posterior_high(timePts3, :)']',1);
    
    % Late
    % % Low
    timePts1 = (X_low == n & (es.outcome(t_low) == 3));
    meanPost.low.late(:,n) = nanmean(Posterior_low(timePts1, :),1);
    % % Norm
    timePts2 = (X_norm == n & (es.outcome(t) == 3));
    meanPost.norm.late(:,n) = nanmean(Posterior_norm(timePts2, :),1);
    % % High
    timePts3 = (X_high == n & (es.outcome(t_high) == 3));
    meanPost.high.late(:,n) = nanmean(Posterior_high(timePts3, :),1);
    % % All
    meanPost.all.late(:,n) = nanmean([Posterior_low(timePts1, :)' Posterior_norm(timePts2, :)' Posterior_high(timePts3, :)']',1);
    
    % Miss
    % % Low
    timePts1 = (X_low == n & (es.outcome(t_low) == 4));
    meanPost.low.miss(:,n) = nanmean(Posterior_low(timePts1, :),1);
    % % Norm
    timePts2 = (X_norm == n & (es.outcome(t) == 4));
    meanPost.norm.miss(:,n) = nanmean(Posterior_norm(timePts2, :),1);
    % % High
    timePts3 = (X_high == n & (es.outcome(t_high) == 4));
    meanPost.high.miss(:,n) = nanmean(Posterior_high(timePts3, :),1);
    % % All
    meanPost.all.miss(:,n) = nanmean([Posterior_low(timePts1, :)' Posterior_norm(timePts2, :)' Posterior_high(timePts3, :)']',1);
    
end

%%
Posterior_all.outcome = outcome;

Posterior_all.meanPost.norm = meanPost_norm;
Posterior_all.meanPost.high = meanPost_high;
Posterior_all.meanPost.low =  meanPost_low;

Posterior_all.meanPost.high_orig = meanPost_high_orig;
Posterior_all.meanPost.low_orig =  meanPost_low_orig;

Posterior_all.meanPost.norm_pass = meanPost_norm_pass;
Posterior_all.meanPost.high_pass = meanPost_high_pass;
Posterior_all.meanPost.low_pass =  meanPost_low_pass;

Posterior_all.meanPost.norm_actv = meanPost_norm_actv;
Posterior_all.meanPost.high_actv = meanPost_high_actv;
Posterior_all.meanPost.low_actv =  meanPost_low_actv;

Posterior_all.meanPost_new.low.early    =  meanPost.low.early;
Posterior_all.meanPost_new.low.correct  =  meanPost.low.correct;
Posterior_all.meanPost_new.low.late     =  meanPost.low.late;
Posterior_all.meanPost_new.low.miss     =  meanPost.low.miss;

Posterior_all.meanPost_new.norm.early    =  meanPost.norm.early;
Posterior_all.meanPost_new.norm.correct  =  meanPost.norm.correct;
Posterior_all.meanPost_new.norm.late     =  meanPost.norm.late;
Posterior_all.meanPost_new.norm.miss     =  meanPost.norm.miss;

Posterior_all.meanPost_new.high.early    =  meanPost.high.early;
Posterior_all.meanPost_new.high.correct  =  meanPost.high.correct;
Posterior_all.meanPost_new.high.late     =  meanPost.high.late;
Posterior_all.meanPost_new.high.miss     =  meanPost.high.miss;

Posterior_all.meanPost_new.all.early    =  meanPost.all.early;
Posterior_all.meanPost_new.all.correct  =  meanPost.all.correct;
Posterior_all.meanPost_new.all.late     =  meanPost.all.late;
Posterior_all.meanPost_new.all.miss     =  meanPost.all.miss;

bestModel = dec.model.bestModel;
bestModel(find(sum(60*bestModel')<1),:) = [];
[~,maxPos] = max(bestModel');
[~,sort_order] = sort(maxPos);
for n = 1:size(bestModel,1)
    bestModel(n,:) = bestModel(n,:) - min(bestModel(n,:));
    bestModel(n,:) = bestModel(n,:) ./max(bestModel(n,:));
end

Posterior_all.data      = es;

if plot_all
    %% Plotting from here
    figure(5)
    subplot(311)
    hold off;
    imagesc(Posterior_low')
    RedWhiteBlue;
    hold on
    axis xy
    plot(1:sum(t_low),X_low,'k.', 'linewidth',2)
    plot(find(es.lick(t_low)==1),X_low(es.lick(t_low)==1),'mo')
    plot(find(es.lick(t_low)==1),X_low(es.lick(t_low)==1),'g.')
    plot(find(es.lick(t_low)==1),X_low(es.lick(t_low)==1),'mo')
    subplot(312)
    hold off;
    imagesc(Posterior_norm')
    axis xy
    hold on;
    plot(1:sum(t),X_norm,'k.', 'linewidth',2)
    plot(find(es.lick(t)==1),X_norm(es.lick(t)==1),'g.')
    plot(find(es.lick(t)==1),X_norm(es.lick(t)==1),'mo')
    subplot(313)
    hold off
    imagesc(Posterior_high')
    axis xy
    hold on;
    plot(1:sum(t_high),X_high,'k.', 'linewidth',2)
    plot(find(es.lick(t_high)==1),X_high(es.lick(t_high)==1),'g.')
    plot(find(es.lick(t_high)==1),X_high(es.lick(t_high)==1),'mo')
    for n = 1:3
        subplot(3,1,n)
        set(gca,'CLim',[-.5 .5])
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    end
    set(gca,'CLim',[-.5 .5])
    hcb = colorbar('YTick',[-.5 0 .5],'YTickLabel',{'2^-0.5 x', 'Chance', '2^0.5 x'});
    set(hcb,'YTickMode','manual')
    % colormap(gray);
    %%
    figure(1)
    imagesc(bestModel(sort_order,:)); axis xy;
    
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    ylabel('Cell #')
    xlabel('Position in room')
    %%
    figure(2)
    subplot(311)
    imagesc(bins_low, dec.bins, meanPost_low); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['Low ' type], 'fontsize', 14);
    subplot(312)
    imagesc(dec.bins, dec.bins, meanPost_norm); axis xy; colorbar; axis tight; axis equal; axis tight
    title('Normal', 'fontsize', 14);
    subplot(313)
    imagesc(bins_high, dec.bins, meanPost_high); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['High ' type], 'fontsize', 14);
    RedWhiteBlue;
    for n = 1:3
        subplot(3,1,n)
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
        line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
    end
    for n = 1:3
        subplot(3,1,n)
        xlabel('Original Position')
        ylabel('Decoded Posterior')
        set(gca,'CLim',[-0.5 0.5])
        hcb = colorbar('YTick',[-1 0 1],'YTickLabel',{'2^-0.5 x', 'Chance', '2^0.5 x'});
        set(hcb,'YTickMode','manual')
    end
    
    %%
    figure(8)
    subplot(321)
    imagesc(dec.bins, dec.bins, meanPost_low_pass); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['Low ' type ': MISS'], 'fontsize', 14);
    RedWhiteBlue;
    subplot(322)
    imagesc(dec.bins, dec.bins, meanPost_low_actv); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['Low ' type ': HIT'], 'fontsize', 14);
    subplot(323)
    imagesc(dec.bins, dec.bins, meanPost_norm_pass); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['Normal: MISS'], 'fontsize', 14);
    subplot(324)
    imagesc(dec.bins, dec.bins, meanPost_norm_actv); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['Normal: HIT'], 'fontsize', 14);
    subplot(325)
    imagesc(dec.bins, dec.bins, meanPost_high_pass); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['High ' type ': MISS'], 'fontsize', 14);
    subplot(326)
    imagesc(dec.bins, dec.bins, meanPost_high_actv); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['High ' type ': HIT'], 'fontsize', 14);
    for n = 1:6
        subplot(3,2,n)
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
        line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
    end
    for n = 1:6
        subplot(3,2,n)
        xlabel('Original Position')
        ylabel('Decoded Posterior')
        set(gca,'CLim',[-0.5 0.5])
        hcb = colorbar('YTick',[-1 0 1],'YTickLabel',{'2^-0.5 x', 'Chance', '2^0.5 x'});
        set(hcb,'YTickMode','manual')
    end
    
    t_n = zeros(size(Posterior_all.Posterior_norm));
    t_l = zeros(size(Posterior_all.Posterior_low));
    t_h = zeros(size(Posterior_all.Posterior_high));
    for n = 1:size(t_l,1)
        t_l(n,Posterior_all.X_low(n)) = 1;
    end
    for n = 1:size(t_h,1)
        t_h(n,Posterior_all.X_high(n)) = 1;
    end
    %     X_norm = X_norm(1:size(Posterior_all.Posterior_norm,1));
    for n = 1:size(t_n,1)
        t_n(n,Posterior_all.X_norm(n)) = 1;
    end
    Posterior_all.confidence.low = (max(Posterior_all.Posterior_low') - min(Posterior_all.Posterior_low'));
    Posterior_all.confidence.norm = (max(Posterior_all.Posterior_norm')- min(Posterior_all.Posterior_norm'));
    Posterior_all.confidence.high = (max(Posterior_all.Posterior_high') - min(Posterior_all.Posterior_high'));
    Posterior_all.confidence.gray = (max(Posterior_all.Posterior_gray') - min(Posterior_all.Posterior_gray'));
    
    Posterior_all.accuracy.low =  Posterior_all.Posterior_low(t_l>0)./max(Posterior_all.Posterior_low(t_l>0)');
    Posterior_all.accuracy.norm = Posterior_all.Posterior_norm(t_n>0)./max(Posterior_all.Posterior_norm(t_n>0)');
    Posterior_all.accuracy.high = Posterior_all.Posterior_high(t_h>0)./max(Posterior_all.Posterior_high(t_h>0)');
    
    [~,X_ML_low] = max(Posterior_all.Posterior_low');
    [~,X_ML_norm] = max(Posterior_all.Posterior_norm');
    [~,X_ML_high] = max(Posterior_all.Posterior_high');
    [~,X_ML_gray] = max(Posterior_all.Posterior_gray');
    
    Posterior_all.MAP.low = X_ML_low;
    Posterior_all.MAP.norm = X_ML_norm;
    Posterior_all.MAP.high = X_ML_high;
    Posterior_all.MAP.gray = X_ML_gray;
    
    Posterior_all.error.low  = ((X_ML_low - X_low'));
    Posterior_all.error.norm = ((X_ML_norm - X_norm'));
    Posterior_all.error.high = ((X_ML_high - X_high'));
    %     Posterior_all.accuracy.low = mean(Posterior_all.Posterior_low(t_l>0));
    %     Posterior_all.accuracy.norm = mean(Posterior_all.Posterior_norm(t>0));
    %     Posterior_all.accuracy.high = mean(Posterior_all.Posterior_high(t_h>0));
    
    %     Posterior_all.width.low = mean(std(Posterior_all.Posterior_low'));
    %     Posterior_all.width.norm = mean(std(Posterior_all.Posterior_norm'));
    %     Posterior_all.width.high = mean(std(Posterior_all.Posterior_high'));
    
    t_l = zeros(size(Posterior_all.Posterior_low_orig));
    t_h = zeros(size(Posterior_all.Posterior_high_orig));
    for n = 1:size(t_l,1)
        t_l(n,Posterior_all.X_low_orig(n)) = 1;
    end
    for n = 1:size(t_h,1)
        t_h(n,Posterior_all.X_high_orig(n)) = 1;
    end
    Posterior_all.confidence.low_orig = (max(Posterior_all.Posterior_low_orig') - min(Posterior_all.Posterior_low_orig'));
    Posterior_all.confidence.high_orig = (max(Posterior_all.Posterior_high_orig') - min(Posterior_all.Posterior_high_orig'));
    
    Posterior_all.accuracy.low_orig = (Posterior_all.Posterior_low_orig(t_l>0))./max(Posterior_all.Posterior_low_orig(t_l>0)');
    Posterior_all.accuracy.high_orig = (Posterior_all.Posterior_high_orig(t_h>0))./max(Posterior_all.Posterior_high_orig(t_h>0)');
    
    [~,X_ML_low_orig] = max(Posterior_all.Posterior_low_orig');
    [~,X_ML_high_orig] = max(Posterior_all.Posterior_high_orig');
    
    Posterior_all.error.low_orig = ((X_ML_low_orig - X_low_orig(1:length(X_ML_low_orig))'));;
    Posterior_all.error.high_orig = ((X_ML_high_orig - X_high_orig(1:length(X_ML_high_orig))'));;
    
    Posterior_all.width.low_orig = mean(std(Posterior_all.Posterior_low_orig'));
    Posterior_all.width.high_orig = mean(std(Posterior_all.Posterior_high_orig'));
    
    %% getting the 45 deg marginals
    tmp = Posterior_all.meanPost.low(:,1:32);
    [Posterior_all.marginals.low, Posterior_all.marginals.lowX] = get45Marginal(tmp);
    tmp = Posterior_all.meanPost.high(:,1:32);
    [Posterior_all.marginals.high, Posterior_all.marginals.highX] = get45Marginal(tmp);
    tmp = Posterior_all.meanPost.norm(:,1:32);
    [Posterior_all.marginals.norm, Posterior_all.marginals.normX] = get45Marginal(tmp);
    
    tmp = Posterior_all.meanPost.low;
    [Posterior_all.marginals.alllow, Posterior_all.marginals.alllowX] = get45Marginal(tmp);
    tmp = Posterior_all.meanPost.high;
    [Posterior_all.marginals.allhigh, Posterior_all.marginals.allhighX] = get45Marginal(tmp);
    tmp = Posterior_all.meanPost.norm;
    [Posterior_all.marginals.allnorm, Posterior_all.marginals.allnormX] = get45Marginal(tmp);
    
    %%
    figure;
    plot(Posterior_all.marginals.lowX*2*sqrt(2),Posterior_all.marginals.low)
    hold on;
    plot(Posterior_all.marginals.normX*2*sqrt(2),Posterior_all.marginals.norm,'k')
    plot(Posterior_all.marginals.highX*2*sqrt(2),Posterior_all.marginals.high,'r');
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    xlabel('Distance from original position (cm)');
    ylabel('Decoded posterior');
    line(xlim, [0 0], 'linestyle','--','color','k');
    line([0 0], ylim, 'linestyle','--','color','k');
    
    %%
    Posterior_all.fieldWidths.norm = 2*(Posterior_all.decoder.numBins - calculate50crosses(Posterior_all.decoder.model.bestModel));
    Posterior_all.fieldWidths.high = 2*(Posterior_all.decoder.numBins - calculate50crosses(Posterior_all.decoder_high.model.bestModel));
    Posterior_all.fieldWidths.low =  2*(Posterior_all.decoder.numBins - calculate50crosses(Posterior_all.decoder_low.model.bestModel));
    
    
    %%
%     figure(15);
%     set(15, 'Position', [100    100   1000   700]);
    f = figure;
    set(f, 'Position', [100    100   1000   700]);
    
    subplot(441)
    imagesc(Posterior_all.meanPost_new.low.early)
    text(3,-5,num2str(length(intersect(outcome.lowContrast, outcome.earlyTrials))));
    subplot(442)
    imagesc(Posterior_all.meanPost_new.low.correct)
    text(3,-5,num2str(length(intersect(outcome.lowContrast, outcome.correctTrials))));
    subplot(443)
    imagesc(Posterior_all.meanPost_new.low.late)
    text(3,-5,num2str(length(intersect(outcome.lowContrast, outcome.lateTrials))));
    subplot(444)
    imagesc(Posterior_all.meanPost_new.low.miss)
    text(3,-5,num2str(length(intersect(outcome.lowContrast, outcome.missTrials))));
    subplot(445)
    imagesc(Posterior_all.meanPost_new.norm.early)
    text(3,-5,num2str(length(intersect(outcome.normContrast, outcome.earlyTrials))));
    subplot(446)
    imagesc(Posterior_all.meanPost_new.norm.correct)
    text(3,-5,num2str(length(intersect(outcome.normContrast, outcome.correctTrials))));
    subplot(447)
    imagesc(Posterior_all.meanPost_new.norm.late)
    text(3,-5,num2str(length(intersect(outcome.normContrast, outcome.lateTrials))));
    subplot(448)
    imagesc(Posterior_all.meanPost_new.norm.miss)
    text(3,-5,num2str(length(intersect(outcome.normContrast, outcome.missTrials))));
    subplot(449)
    imagesc(Posterior_all.meanPost_new.high.early)
    text(3,-5,num2str(length(intersect(outcome.highContrast, outcome.earlyTrials))));
    subplot(4,4,10)
    imagesc(Posterior_all.meanPost_new.high.correct)
    text(3,-5,num2str(length(intersect(outcome.highContrast, outcome.correctTrials))));
    subplot(4,4,11)
    imagesc(Posterior_all.meanPost_new.high.late)
    text(3,-5,num2str(length(intersect(outcome.highContrast, outcome.lateTrials))));
    subplot(4,4,12)
    imagesc(Posterior_all.meanPost_new.high.miss)
    text(3,-5,num2str(length(intersect(outcome.highContrast, outcome.missTrials))));
    
    subplot(4,4,13)
    imagesc(Posterior_all.meanPost_new.all.early)
    text(3,-5,num2str(length(intersect(outcome.allContrast, outcome.earlyTrials))));
    subplot(4,4,14)
    imagesc(Posterior_all.meanPost_new.all.correct)
    text(3,-5,num2str(length(intersect(outcome.allContrast, outcome.correctTrials))));
    subplot(4,4,15)
    imagesc(Posterior_all.meanPost_new.all.late)
    text(3,-5,num2str(length(intersect(outcome.allContrast, outcome.lateTrials))));
    subplot(4,4,16)
    imagesc(Posterior_all.meanPost_new.all.miss)
    text(3,-5,num2str(length(intersect(outcome.allContrast, outcome.missTrials))));
    for n = 1:16
        subplot(4,4,n)
        axis xy; axis equal; axis tight;
        line([1 50], [1 50],'linestyle','--','color','k')
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        set(gca,'CLim',[-0.5 0.5])
        set(gca, 'XTick',[0 25 50],'XTickLabel',[],'YTick',[0 25 50],'YTickLabel',[]);
    end
    subplot(4,4,1)
    ylabel('Low Contrast')
    title('EARLY')
    subplot(4,4,2)
    title('CORRECT')
    subplot(4,4,3)
    title('LATE')
    subplot(4,4,4)
    title('MISS')
    subplot(4,4,5)
    ylabel('Baseline Contrast');
    subplot(4,4,9)
    ylabel('High Contrast');
    subplot(4,4,13)
    ylabel('All');
    RedWhiteBlue;
    drawnow; 
    pause(1);
end