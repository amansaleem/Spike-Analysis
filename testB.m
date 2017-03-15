for icell = 1:length(es.spikeIDs)
    for ibin = 1:10
        bins = (es.traj>((ibin-1)*10) & es.ballspeed>1 & es.traj<=(ibin*10));
        FR(icell,ibin) = nanmean(es.spikeTrain(bins,icell))*60;
    end
    for ishuf = 1:25
        s = es.spikeTrain(randperm(size(es.spikeTrain,1)),icell);
        for ibin = 1:10
            bins = (es.traj>((ibin-1)*10) & es.ballspeed>1 & es.traj<=(ibin*10));
            FRs(ishuf).FR(icell,ibin) = nanmean(s(bins))*60;
        end
    end
end

figure;
idx = 1;
for icell = [1 2 3 5 6 9 11 12 13 17 23 24 26 27 28 30]
    subplot(4,4,idx)
    idx = idx + 1;
    hold on;
    for ishuf = 1:25
        plot(5:10:95,FRs(ishuf).FR(icell,:),'color',[0.5 0.5 0.5]);
    end
    plot(5:10:95,FR(icell,:),'.-r','linewidth',1.5);
    axis tight;
    
    line([20 20],ylim,'color',[.5 .5 .5],'linestyle','--')
    line([40 40],ylim,'color',[.5 .5 .5],'linestyle','--')
    line([60 60],ylim,'color',[.5 .5 .5],'linestyle','--')
    line([80 80],ylim,'color',[.5 .5 .5],'linestyle','--')
    
    title(['Cell: ' num2str(icell)]);
    set(gca, 'TickDir','out', 'box','off','color','none','fontsize',14,'XTick',[]);
end

figure;
for icell = 28
    plot(es.sampleTimes, es.traj,'k')
    hold on;
    k = find(es.spikeTrain(1:length(es.traj),icell) & es.ballspeed>2);
    plot(es.sampleTimes(k), es.traj(k),'r.');
    title(num2str(icell));
    pause
    hold off;
end