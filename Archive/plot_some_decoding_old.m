clear all
% load con_920_2
% es = con_920.data.es;
load con_918_2
es = con_918.data.es;

% load con_918_V1
% es = con_918_V1.data.es;
%%
test_normal = 1;

base = es.traj~=0 & ~isnan(es.traj) & es.gain==1 & es.roomLength==1 & es.smthBallSpd>5;

% base = base & (es.iexp<104);

t      = es.contrast==0.6 & base;
t_low  = es.contrast==0.18 & base;
t_high = es.contrast==0.72 & base;

dec = bayesDecoder;
dec.numBins = 50;
smthWin = 200;

%%
for icell = 1:size(es.spikeTrain,2);
    spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), 60, smthWin);
    
    %     spkRate_norm = smthInTime(es.spikeTrain(:,icell),60,500, 'same', t);
    %     spkRate_low  = smthInTime(es.spikeTrain(:,icell),60,500, 'same', t_low);
    %     spkRate_high = smthInTime(es.spikeTrain(:,icell),60,500, 'same', t_high);
    %
    %     spkRate(t     ,icell)     = spkRate_norm(t);
    %     spkRate(t_low ,icell)     = spkRate_norm(t_low);
    %     spkRate(t_high,icell)     = spkRate_norm(t_high);
    
end
if ~test_normal
    dec_l = dec;
    dec_h = dec;
end
[dec, ~, X_norm, Posterior_norm, nPosterior_norm] = ...
    dec.trainDecoder( es.traj(t), spkRate(t,:), 0);

if ~test_normal
    [dec_l, ~, X_low, Posterior_low, nPosterior_low] = ...
        dec_l.trainDecoder( es.traj(t_low), spkRate(t_low,:), 0);
    [dec_h, ~, X_high, Posterior_high, nPosterior_high] = ...
        dec_h.trainDecoder( es.traj(t_high), spkRate(t_high,:), 0);
else
    [~,  Posterior_low, nPosterior_low] = dec.predictBayesDecoder(spkRate(t_low,:), 0,'best');
    [~,  Posterior_high, nPosterior_high] = dec.predictBayesDecoder(spkRate(t_high,:), 0,'best');
    % [~,  Posterior_norm, nPosterior_norm] = dec.predictBayesDecoder(spkRate(t,:), 0,'best');
end

Posterior_low   = log2(Posterior_low*dec.numBins);
Posterior_high  = log2(Posterior_high*dec.numBins);
nPosterior_low   = log2(nPosterior_low*dec.numBins);
nPosterior_high  = log2(nPosterior_high*dec.numBins);

Posterior_norm  = log2(Posterior_norm*dec.numBins);
nPosterior_norm  = log2(nPosterior_norm*dec.numBins);

if test_normal
    X_low = normalise1var(es.traj(t_low), dec.numBins);
    % X_norm = normalise1var(es.traj(t), 50);
    X_high = normalise1var(es.traj(t_high), dec.numBins);
else
    X_high = X_high(1:size(Posterior_high,1));
    X_low = X_low(1:size(Posterior_low,1));
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
%%
for n = 1:dec.numBins
    timePts = (X_norm == n);
    meanPost_norm(:,n) = nanmean(Posterior_norm(timePts, :),1);
    timePts = (X_high == n);
    meanPost_high(:,n) = nanmean(Posterior_high(timePts, :),1);
    timePts = (X_low == n);
    meanPost_low(:,n) = nanmean(Posterior_low(timePts, :),1);
    timePts = (X_norm == n & (es.outcome(t) ~= 2));
    meanPost_norm_pass(:,n) = nanmean(Posterior_norm(timePts, :),1);
    timePts = (X_norm == n & (es.outcome(t) == 2));
    meanPost_norm_actv(:,n) = nanmean(Posterior_norm(timePts, :),1);
    timePts = (X_low == n & (es.outcome(t_low) == 2));
    meanPost_low_actv(:,n) = nanmean(Posterior_low(timePts, :),1);
    timePts = (X_low == n & (es.outcome(t_low) ~= 2));
    meanPost_low_pass(:,n) = nanmean(Posterior_low(timePts, :),1);
    timePts = (X_high == n & (es.outcome(t_high) ~= 2));
    meanPost_high_pass(:,n) = nanmean(Posterior_high(timePts, :),1);
    timePts = (X_high == n & (es.outcome(t_high) == 2));
    meanPost_high_actv(:,n) = nanmean(Posterior_high(timePts, :),1);
end



bestModel = dec.model.bestModel;
bestModel(find(sum(60*bestModel')<1),:) = [];
[~,maxPos] = max(bestModel');
[~,sort_order] = sort(maxPos);
for n = 1:size(bestModel,1)
    bestModel(n,:) = bestModel(n,:) - min(bestModel(n,:));
    bestModel(n,:) = bestModel(n,:) ./max(bestModel(n,:));
end

%%
figure(5)
subplot(311)
hold off;
imagesc(Posterior_low')
hold on
axis xy
plot(1:sum(t_low),X_low,'k')
plot(find(es.lick(t_low)==1),X_low(es.lick(t_low)==1),'mo')
plot(find(es.lick(t_low)==1),X_low(es.lick(t_low)==1),'g.')
plot(find(es.lick(t_low)==1),X_low(es.lick(t_low)==1),'mo')
subplot(312)
hold off;
imagesc(Posterior_norm')
axis xy
hold on;
plot(1:sum(t),X_norm,'k')
plot(find(es.lick(t)==1),X_norm(es.lick(t)==1),'g.')
plot(find(es.lick(t)==1),X_norm(es.lick(t)==1),'mo')
subplot(313)
hold off
imagesc(Posterior_high')
axis xy
hold on;
plot(1:sum(t_high),X_high,'k')
plot(find(es.lick(t_high)==1),X_high(es.lick(t_high)==1),'g.')
plot(find(es.lick(t_high)==1),X_high(es.lick(t_high)==1),'mo')
for n = 1:3
    subplot(3,1,n)
    set(gca,'CLim',[-1 1])
end
set(gca,'CLim',[-1 1])
hcb = colorbar('YTick',[-1 0 1],'YTickLabel',{'1/2 x', 'Chance', '2x'});
set(hcb,'YTickMode','manual')
%%
figure(1)
imagesc(bestModel(sort_order,:)); axis xy;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
ylabel('Cell #')
xlabel('Position in room')
%%
figure(2)
subplot(311)
imagesc(dec.bins, dec.bins, meanPost_low_pass); axis xy; colorbar; axis tight; axis equal; axis tight
title('Low contrast', 'fontsize', 14);
subplot(312)
imagesc(dec.bins, dec.bins, meanPost_norm); axis xy; colorbar; axis tight; axis equal; axis tight
title('Normal', 'fontsize', 14);
subplot(313)
imagesc(dec.bins, dec.bins, meanPost_high); axis xy; colorbar; axis tight; axis equal; axis tight
title('High contrast', 'fontsize', 14);
for n = 1:3
    subplot(3,1,n)
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
end
for n = 1:3
    subplot(3,1,n)
    xlabel('Original Position')
    ylabel('Decoded Posterior')
    set(gca,'CLim',[-1 1])
    hcb = colorbar('YTick',[-1 0 1],'YTickLabel',{'1/2 x', 'Chance', '2x'});
    set(hcb,'YTickMode','manual')
end

%%
figure(8)
subplot(321)
imagesc(dec.bins, dec.bins, meanPost_low_pass); axis xy; colorbar; axis tight; axis equal; axis tight
title('Low contrast: MISS', 'fontsize', 14);
subplot(322)
imagesc(dec.bins, dec.bins, meanPost_low_actv); axis xy; colorbar; axis tight; axis equal; axis tight
title('Low contrast: HIT', 'fontsize', 14);
subplot(323)
imagesc(dec.bins, dec.bins, meanPost_norm_pass); axis xy; colorbar; axis tight; axis equal; axis tight
title('Normal: MISS', 'fontsize', 14);
subplot(324)
imagesc(dec.bins, dec.bins, meanPost_norm_actv); axis xy; colorbar; axis tight; axis equal; axis tight
title('Normal: HIT', 'fontsize', 14);
subplot(325)
imagesc(dec.bins, dec.bins, meanPost_high_pass); axis xy; colorbar; axis tight; axis equal; axis tight
title('High contrast: MISS', 'fontsize', 14);
subplot(326)
imagesc(dec.bins, dec.bins, meanPost_high_actv); axis xy; colorbar; axis tight; axis equal; axis tight
title('High contrast: HIT', 'fontsize', 14);
for n = 1:6
    subplot(3,2,n)
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
end
for n = 1:6
    subplot(3,2,n)
    xlabel('Original Position')
    ylabel('Decoded Posterior')
    set(gca,'CLim',[-1 1])
    hcb = colorbar('YTick',[-1 0 1],'YTickLabel',{'1/2 x', 'Chance', '2x'});
    set(hcb,'YTickMode','manual')
end

%%
figure(10)

subplot(311)
imagesc(dec.bins, dec.bins, meanPost_low_actv - meanPost_low_pass); axis xy; colorbar; axis tight; axis equal; axis tight
title('Low contrast: Hit - Miss', 'fontsize', 14);
subplot(312)
imagesc(dec.bins, dec.bins, meanPost_norm_actv - meanPost_norm_pass); axis xy; colorbar; axis tight; axis equal; axis tight
title('Normal: Hit - Miss', 'fontsize', 14);
subplot(313)
imagesc(dec.bins, dec.bins, meanPost_high_actv - meanPost_high_pass); axis xy; colorbar; axis tight; axis equal; axis tight
title('High contrast: Hit - Miss', 'fontsize', 14);

for n = 1:3
    subplot(3,1,n)
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
end
for n = 1:3
    subplot(3,1,n)
    xlabel('Original Position')
    ylabel('Decoded Posterior')
    set(gca,'CLim',[-1 1])
    hcb = colorbar('YTick',[-1 0 1],'YTickLabel',{'1/2 x', 'Chance', '2x'});
    set(hcb,'YTickMode','manual')
end

%%
figure(11)
bestModel = dec.model.bestModel;
[~,maxPos] = max(bestModel');
reward_cell = (maxPos<=35 & maxPos>=28);

for icell = 1:size(es.spikeTrain,2)
    hit_rewpos_low(icell) = nanmean(spkRate(es.outcome(t_low)==2 & X_low>=28 & X_low<=35, icell));
    mis_rewpos_low(icell) = nanmean(spkRate(es.outcome(t_low)~=2 & X_low>=28 & X_low<=35, icell));
    
    hit_rewpos_high(icell) = nanmean(spkRate(es.outcome(t_high)==2 & X_high>=28 & X_high<=35, icell));
    mis_rewpos_high(icell) = nanmean(spkRate(es.outcome(t_high)~=2 & X_high>=28 & X_high<=35, icell));
    
    hit_rewpos_norm(icell) = nanmean(spkRate(es.outcome(t)==2 & X_norm>=28 & X_norm<=35, icell));
    mis_rewpos_norm(icell) = nanmean(spkRate(es.outcome(t)~=2 & X_norm>=28 & X_norm<=35, icell));
    
    hit_Nrewpos_low(icell) = nanmean(spkRate(es.outcome(t_low)==2 & ~(X_low>=28 & X_low<=35), icell));
    mis_Nrewpos_low(icell) = nanmean(spkRate(es.outcome(t_low)~=2 & ~(X_low>=28 & X_low<=35), icell));
    
    hit_Nrewpos_high(icell) = nanmean(spkRate(es.outcome(t_high)==2 & ~(X_high>=28 & X_high<=35), icell));
    mis_Nrewpos_high(icell) = nanmean(spkRate(es.outcome(t_high)~=2 & ~(X_high>=28 & X_high<=35), icell));
    
    mis_Nrewpos_norm(icell) = nanmean(spkRate(es.outcome(t)~=2 & ~(X_norm>=28 & X_norm<=35), icell));
    hit_Nrewpos_norm(icell) = nanmean(spkRate(es.outcome(t)==2 & ~(X_norm>=28 & X_norm<=35), icell));
end

subplot(321)
plot(60*mis_rewpos_low, 60*hit_rewpos_low, 'k.'); axis equal;
hold on;
plot(60*mis_rewpos_low(reward_cell), 60*hit_rewpos_low(reward_cell), 'r.'); axis equal;
axis([0 40 0 40])
line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
xlabel('Rate on hits')
xlabel('Rate on misses')
ylabel('Rate on Hits')

subplot(323)
plot(60*mis_rewpos_norm, 60*hit_rewpos_norm, 'k.'); axis equal;
hold on;
plot(60*mis_rewpos_norm(reward_cell), 60*hit_rewpos_norm(reward_cell), 'r.'); axis equal;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis([0 40 0 40])
line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')

subplot(325)
plot(60*mis_rewpos_high, 60*hit_rewpos_high, 'k.'); axis equal;
hold on;
plot(60*mis_rewpos_high(reward_cell), 60*hit_rewpos_high(reward_cell), 'r.'); axis equal;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis([0 40 0 40])
line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')

subplot(322)
plot(60*hit_rewpos_low, 60*hit_Nrewpos_low, 'k.'); axis equal;
hold on;
plot(60*hit_rewpos_low(reward_cell), 60*hit_Nrewpos_low(reward_cell), 'r.'); axis equal;
axis([0 40 0 40])
line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
xlabel('Rate on hits')
xlabel('Rate on misses')
ylabel('Rate on Hits')

subplot(324)
plot(60*hit_rewpos_norm, 60*hit_Nrewpos_norm, 'k.'); axis equal;
hold on;
plot(60*hit_rewpos_norm(reward_cell), 60*hit_Nrewpos_norm(reward_cell), 'r.'); axis equal;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis([0 40 0 40])
line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')

subplot(326)
plot(60*hit_rewpos_high, 60*hit_Nrewpos_high, 'k.'); axis equal;
hold on;
plot(60*hit_rewpos_high(reward_cell), 60*hit_Nrewpos_high(reward_cell), 'r.'); axis equal;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
axis([0 40 0 40])
line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
