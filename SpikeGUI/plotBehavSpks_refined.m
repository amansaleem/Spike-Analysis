function plotBehavSpks_refined(es,t,var,plot_index,icell,delayT)

if nargin<5
    showSpikes = 0;
else
    showSpikes = 1;
end
if nargin<6
    delayT = 0;
end
if showSpikes
    
    es.spikeTrain = circshift(es.spikeTrain,[-delayT 0]);
    
    % Taking out sections of the data where there are no spikes at all
    if sum(es.spikeTrain(:,icell))>20
        spkTrials_start = 1;
        spkTrials_end   = max(es.trialID);
        while sum(es.spikeTrain(es.trialID==spkTrials_start,icell))==0 | spkTrials_start==spkTrials_end
            spkTrials_start = spkTrials_start + 1;
        end
        while sum(es.spikeTrain(es.trialID==spkTrials_end,icell))==0  | spkTrials_start==spkTrials_end
            spkTrials_end = spkTrials_end - 1;
        end
        goodTrials = ones(size(es.traj));
        goodTrials(es.trialID < spkTrials_start) = 0;
        goodTrials(es.trialID > spkTrials_end) = 0;
        
        t = t & goodTrials>0;
    end
end
%
if size(var,2)==1
    trialIDs = unique(es.trialID(t));
    con = nanmean(es.contrast(t));
    tmp_trialID = nan*ones(size(es.trialID));
    for itr = 1:length(trialIDs)
        %% make shaded regions
        ROI = es.trialID==trialIDs(itr);
        tmp_trialID(ROI) = itr;
        minTraj = min(es.traj(ROI));
        maxTraj = max(es.traj(ROI));
        X = [minTraj maxTraj maxTraj minTraj];
        Y = [itr-1   itr-1   itr     itr    ];
        switch nanmean(es.outcome(ROI))
            case 0 % Passive, light gray
                patch(X,Y,[0.1 0.1 0.1], ...
                    'EdgeColor','none','linestyle','none');
            case 1 % Early, light red
                patch(X,Y,[1 0.65 0.65], ...
                    'EdgeColor','none','linestyle','none');
            case 2 % Correct, gray
                patch(X,Y,[0.5 0.5 0.5], ...
                    'EdgeColor','none','linestyle','none');
            case 3 % Late, light Blue
                patch(X,Y,[0.65 0.65 1], ...
                    'EdgeColor','none','linestyle','none');
            case 4 % Miss, white
                patch(X,Y,[1 1 1], ...
                    'EdgeColor','none','linestyle','none');
            case 5 % Time-out, light gray
                patch(X,Y,[0.1 0.1 0.1], ...
                    'EdgeColor','none','linestyle','none');
        end
        %         plot(plot_index,var(tmp_trialID==itr),tmp_trialID(tmp_trialID==itr), 'color',[.6 .6 .6])
        hold(plot_index,'on');
    end
    es.trialID = tmp_trialID;
    % Need to move this
    %     plot(plot_index,var(t),es.trialID(t), 'color',[.7 .7 .7])
    
    if showSpikes
        if ~isfield(es,'theta')
            %             plot(plot_index,es.trialID(es.lick & t), var(es.lick & t),'.','color',[.7 .7 .7])
            plot(plot_index, var(es.spikeTrain(:,icell)>0 & t & es.outcome~=2), es.trialID(es.spikeTrain(:,icell)>0 & t & es.outcome~=2),'rs')
            plot(plot_index, var(es.spikeTrain(:,icell)>0 & t & es.outcome==2), es.trialID(es.spikeTrain(:,icell)>0 & t & es.outcome==2),'bs')
        else
            hold(plot_index,'off');
            colormap_val = colormap(hsv);
            %             plot(plot_index,es.trialID(es.lick & t), var(es.lick & t),'.','color',[.7 .7 .7])
            %             for ispike = find(es.spikeTrain(:,icell) & t)'
            %                 plot(plot_index,es.trialID(ispike), var(ispike),'Marker','+','color',colormap_val(32+round(62*(phase(es.theta.A.hill(ispike)))./(2*pi)),:));
            plot(plot_index, var(es.spikeTrain(:,icell) & t & es.outcome==2), ...
                180+(phase(es.theta.B.hill(es.spikeTrain(:,icell) & t & es.outcome==2))).*(360/(2*pi)), 'k.','MarkerSize',5);
            hold(plot_index,'on');
            %             plot(plot_index, var(es.spikeTrain(:,icell) & t & es.outcome==2), ...
            %                 180+360+(phase(es.theta.B.hill(es.spikeTrain(:,icell) & t & es.outcome==2))).*(360/(2*pi)), 'k.');
            plot(plot_index, var(es.spikeTrain(:,icell) & t & es.outcome==2), ...
                180+360+(phase(es.theta.B.hill(es.spikeTrain(:,icell) & t & es.outcome==2))).*(360/(2*pi)), '.','color',-.5+[.5 .5 .5],'MarkerSize',5);
            %             plot(plot_index, var(es.spikeTrain(:,icell) & t & es.outcome==2), ...
            %                 360+40+360+180+360+(phase(es.theta.A.hill(es.spikeTrain(:,icell) & t & es.outcome==2))).*(360/(2*pi)), '.','color',[.5 .5 .5]);
            %             end
            set(plot_index, 'XLim', [min(var(t)) max(var(t))]);
            set(plot_index, 'YLim', [0 360*2]);
        end
    else
        if ~isfield(es,'theta')
            axis tight;
            %% only licks, rewards & shaded regions
            line([es.rewardPos(10)-es.rewardTolerance es.rewardPos(10)-es.rewardTolerance],...
                ylim, 'color','c');
            line([es.rewardPos(10)+es.rewardTolerance es.rewardPos(10)+es.rewardTolerance],...
                ylim, 'color','c');
            likTimes = es.lick & t;
            rewTimes = es.reward>0 & t;
            line([var(likTimes) var(likTimes)]', ...
                [es.trialID(likTimes)-1 es.trialID(likTimes)]',...
                'linewidth',1,'color','k');
            plot(plot_index, var(rewTimes), es.trialID(rewTimes)-0.5,'co', 'MarkerFaceColor', 'c')
        else
            show_color = 0;
            if show_color>0
                plot(plot_index, var(es.lick & t & es.outcome==2),...
                    180+(phase(es.theta.B.hill(es.lick & t & es.outcome==2))).*(360/(2*pi)),'r+', 'MarkerSize',3)
                hold(plot_index,'on');
                plot(plot_index, var(es.lick & t & es.outcome==2), ...
                    360+180+(phase(es.theta.B.hill(es.lick & t & es.outcome==2))).*(360/(2*pi)),'r+', 'MarkerSize',3)
                plot(plot_index,es.trialID(es.lick & t), var(es.lick & t),'.','color',[.7 .7 .7])
            else
                hold(plot_index,'off');
                colormap_val = colormap(hsv);
                for ilick = find(es.lick & t & es.outcome==2)'
                    plot(plot_index,var(ilick),es.trialID(ilick), 'Marker','+','color',colormap_val(32+round(62*(phase(es.theta.A.hill(ilick)))./(2*pi)),:));
                    hold(plot_index,'on');
                end
            end
        end
    end
    
    set(plot_index, 'box','off','TickDir','out','fontsize',14,'color','none')
    hold(plot_index,'off');
    axis(plot_index,'tight');
else
    plot(plot_index,var(t,2), var(t,1),'color',[.75 .75 .75])
    hold(plot_index,'on');
    if showSpikes
        %         plot(plot_index,var(es.lick & t,2), var(es.lick & t,1),'.','color',[.7 .7 .7])
        
        if ~isfield(es,'A')
            %             plot(plot_index,es.trialID(es.lick & t), var(es.lick & t),'.','color',[.7 .7 .7])
            plot(plot_index,var(es.spikeTrain(:,icell)>0 & t & es.outcome~=2,2), var(es.spikeTrain(:,icell)>0 & t & es.outcome~=2,1),'rs');
            plot(plot_index,var(es.spikeTrain(:,icell)>0 & t & es.outcome==2,2), var(es.spikeTrain(:,icell)>0 & t & es.outcome==2,1),'bs');
        else
            colormap_val = colormap(hsv);
            %             plot(plot_index,es.trialID(es.lick & t), var(es.lick & t),'.','color',[.7 .7 .7])
            for ispike = find(es.spikeTrain(:,icell) & t)'
                plot(plot_index,var(ispike,2), var(ispike,1),'Marker','+','color',colormap_val(32+round(62*(phase(es.theta.A.hill(ispike)))./(2*pi)),:));
            end
        end
    else
        plot(plot_index,var(es.lick & t,2), var(es.lick & t,1),'r+')
        plot(plot_index,var(es.reward>0 & t,2), var(es.reward>0 & t,1),'ko')
    end
    
    set(plot_index, 'box','off','TickDir','out','fontsize',14,'color','none')
    hold(plot_index,'off');
    axis(plot_index,'tight');
end
end