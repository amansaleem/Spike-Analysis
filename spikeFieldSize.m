% function run_all_decoding
run_series = 1:8;
idx = 0;
smthWin = 250;
box_filt = 1;
trainCorrect = 0;
quickProcess = 0;


idx = idx + 1;
expt_list(idx).animal     = 'M130920_BALL';
expt_list(idx).iseries    = 1025;
expt_list(idx).expt_list  = 102:103;
expt_list(idx).shank_list  = 0:7;

idx = idx + 1;
expt_list(idx).animal     = 'M130918_BALL';
expt_list(idx).iseries    = 1030;
expt_list(idx).expt_list  = 103:105;
expt_list(idx).shank_list  = 0:7;

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 530;
expt_list(idx).expt_list  = 104:106;
expt_list(idx).shank_list  = 0:7;

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 531;
expt_list(idx).expt_list  = 103:106;
expt_list(idx).shank_list  = 0:7;

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 601;
expt_list(idx).expt_list  = 103:106;
expt_list(idx).shank_list  = [0:3 5:7];

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 602;
expt_list(idx).expt_list  = 102:106;
expt_list(idx).shank_list  = 0:7;

idx = idx + 1;
expt_list(idx).animal     = 'M140502_BALL';
expt_list(idx).iseries    = 603;
expt_list(idx).expt_list  = 107:110;
expt_list(idx).shank_list  = 2:7;

idx = idx + 1;
expt_list(idx).animal     = 'M140502_BALL';
expt_list(idx).iseries    = 604;
expt_list(idx).expt_list  = 107:110;
expt_list(idx).shank_list  = [1 3 5 6];

% skagg_n = [];
% EV_n = [];
% mRate_n = [];
% skagg_h = [];
% EV_h = [];
% mRate_h = [];
% skagg_l = [];
% EV_l = [];
% mRate_l = [];
% zskagg_n_all = [];
% zskagg_l_all = [];
% zskagg_h_all = [];

for idx = run_series
    expt_list(idx).iseries
    SetDirs;
    es = VRLoadMultipleExpts(expt_list(idx).animal, expt_list(idx).iseries, expt_list(idx).expt_list,'SPIKES',[],expt_list(idx).shank_list);
    
    switch expt_list(idx).iseries
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
    
    base = es.traj~=0 & es.contrast~=0 & ~isnan(es.traj) ...
        & es.gain==1 & es.roomLength==1 & es.smthBallSpd>5 ...
        & es.trajspeed>=0; % & es.traj<=65;
    t      = es.contrast==0.6 & base;
    t_low  = es.contrast==0.18 & base;
    t_high = es.contrast==0.72 & base;
    t_gray = es.traj~=0 & es.contrast==0 & ~isnan(es.traj) & es.smthBallSpd>5;
    
    spkRate = zeros(size(es.spikeTrain));
    es.spikeTrain(es.traj==0,:) = nan; % to see if there is an effect of the gray screen.
    
    for icell = 1:size(es.spikeTrain,2);
        if box_filt
            spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), 60, smthWin,'box');
        else
            spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), 60, smthWin); %*(smthWin*60./1000);
        end
    end
    
    map = oneDimMap;
    [map, pred, X] = map.trainSpikeMap(es.traj(t), spkRate(t,:),0);
    input_data = es.traj(t);
    input_spike = es.spikeTrain(t,:);
        
    map_l = oneDimMap;
    [map_l, pred_l, X_l] = map_l.trainSpikeMap(es.traj(t_low), ...
        spkRate(t_low,:),0);
    input_low = es.traj(t_low);
    input_spike_low = es.spikeTrain(t_low,:);
    
    map_h = oneDimMap;
    [map_h, pred_h, X_h] = map_h.trainSpikeMap(es.traj(t_high), ...
        spkRate(t_high,:),0);
    input_high = es.traj(t_high);
    input_spike_high = es.spikeTrain(t_high,:);
    
    numShuff = 20;
    for iShuff = 1:numShuff
        disp(num2str(iShuff));
        shuffRate_n = input_spike(randperm(length(input_data)),:);
        shuffRate_h = input_spike_high(randperm(length(input_high)),:);
        shuffRate_l = input_spike_low(randperm(length(input_low)),:);
        
        for icell = 1:size(es.spikeTrain,2);
            shuffRate_n(:,icell) = smthInTime(shuffRate_n(:,icell), 60, smthWin,'box');
            shuffRate_h(:,icell) = smthInTime(shuffRate_h(:,icell), 60, smthWin,'box');
            shuffRate_l(:,icell) = smthInTime(shuffRate_l(:,icell), 60, smthWin,'box');
        end
        map_s(iShuff).map = oneDimMap;
        [map_s(iShuff).map] = map_s(iShuff).map.trainSpikeMap( ...
            input_data, shuffRate_n,0);
        
        map_sl(iShuff).map = oneDimMap;
        [map_sl(iShuff).map] = map_sl(iShuff).map.trainSpikeMap( ...
            input_low, shuffRate_l,0);
        
        map_sh(iShuff).map = oneDimMap;
        [map_sh(iShuff).map] = map_sh(iShuff).map.trainSpikeMap( ...
            input_high, shuffRate_h,0);
    end
    
    clear distr_n distr_l distr_h
    for n = 1:length(es.spikeIDs)
        for k = 1:5
            for m = 1:numShuff
                distr_n(n,m,k) = map_s(m).map.model.skaggs(k,n);
                distr_l(n,m,k) = map_sl(m).map.model.skaggs(k,n);
                distr_h(n,m,k) = map_sh(m).map.model.skaggs(k,n);
            end
        end
    end
    actual_n = nanmean(map.model.skaggs,1);
    actual_l = nanmean(map_l.model.skaggs,1);
    actual_h = nanmean(map_h.model.skaggs,1);
    
    distr_n = reshape(distr_n, [size(distr_n,1) size(distr_n,2)*size(distr_n,3)]);
    distr_l = reshape(distr_l, [size(distr_l,1) size(distr_l,2)*size(distr_l,3)]);
    distr_h = reshape(distr_h, [size(distr_h,1) size(distr_h,2)*size(distr_h,3)]);
    
    zskagg_n = (actual_n - nanmean(distr_n,2)')./nanstd(distr_n');
    zskagg_l = (actual_l - nanmean(distr_l,2)')./nanstd(distr_l');
    zskagg_h = (actual_h - nanmean(distr_h,2)')./nanstd(distr_h');
    
    thres = 0.01;
    
    disp(['Norm over threshold: ' num2str(sum(nanmean(map.model.EV)>=thres))])
    
    disp(['All over threshold: ' num2str(sum(nanmean(map.model.EV)>=thres &...
                nanmean(map_h.model.EV)>=thres...
                & nanmean(map_l.model.EV)>=thres))])
        %%   
        for n = 1:size(map.model.EV,2)
        if nanmean(map.model.EV(:,n))>=thres &&...
                nanmean(map_h.model.EV(:,n))>=thres...
                && nanmean(map_l.model.EV(:,n))>=thres
            placeField_norm = 60*nanmean(map.model.tuning(n).respModel,1);
            placeField_high = 60*nanmean(map_h.model.tuning(n).respModel,1);
            placeField_low  = 60*nanmean(map_l.model.tuning(n).respModel,1);
            subplot(121)
            plot(placeField_norm,'k');
            hold on;
            plot(placeField_low,'b');
            plot(placeField_high,'r');
            hold off;
            title([num2str(n) '    NL:'...
                num2str(corr(placeField_norm',placeField_low'))...
                '    NH:' num2str(corr(placeField_norm',placeField_high'))...
                '    HL:' num2str(corr(placeField_high',placeField_low'))])
            legend(['N:' num2str(round(1000*nanmean(map.model.EV(:,n)))/10)],...
                ['L:' num2str(round(1000*nanmean(map_l.model.EV(:,n)))/10)],...
                ['H:' num2str(round(1000*nanmean(map_h.model.EV(:,n)))/10)])
            
            subplot(243)
            plot(zskagg_l,...
                zskagg_h,'ko')
            hold on;
            plot(zskagg_l(n),...
                zskagg_h(n),'r*');
            hold off
            axis equal; axis([0 30 0 30]);
            line(xlim, ylim);
            
            zskagg_n_all = [zskagg_n_all zskagg_n(n)];
            zskagg_h_all = [zskagg_h_all zskagg_h(n)];
            zskagg_l_all = [zskagg_l_all zskagg_l(n)];
            
            skagg_n = [skagg_n actual_n(n)];
            skagg_h = [skagg_h actual_h(n)];
            skagg_l = [skagg_l actual_l(n)];
            
%             hold on;            
%             subplot(243)
            placeField_norm = normalise1var(placeField_norm);
            placeField_high = normalise1var(placeField_high);
            placeField_low  = normalise1var(placeField_low);
%             plot(placeField_norm,'k');
%             hold on;
%             plot(placeField_low,'b');
%             plot(placeField_high,'r');
%             hold off;
            
            subplot(244)
            plot(nanmean(map_l.model.EV,1),...
                nanmean(map_h.model.EV,1),'ko')
            hold on;
            plot(nanmean(map_l.model.EV(:,n)),...
                nanmean(map_h.model.EV(:,n)),'r*')
            hold off;
            axis equal; axis([-0.05 0.3 -0.05 0.3])
            line(xlim, ylim);
            EV_n = [EV_n nanmean(map.model.EV(:,n))];
            EV_h = [EV_h nanmean(map_h.model.EV(:,n))];
            EV_l = [EV_l nanmean(map_l.model.EV(:,n))];
                        
            subplot(247)
            placeField_norm = sort(placeField_norm);
            placeField_high = sort(placeField_high);
            placeField_low  = sort(placeField_low);
            plot(placeField_norm,'k');
            hold on;
            plot(placeField_low,'b');
            plot(placeField_high,'r');
            line(xlim, [0.5 0.5], 'linestyle','--','color','k')
            line([25 25], ylim, 'linestyle','--','color','k')
            legend(['N:' num2str(60*nanmean(spkRate(t,n)))],...
                ['L:' num2str(60*nanmean(spkRate(t_low,n)))],...
                ['H:' num2str(60*nanmean(spkRate(t_high,n)))],...
                'Location','NorthWest')
            hold off;
            
            subplot(248)
            hold on;
            plot(60*nanmean(spkRate(t_low,n)),60*nanmean(spkRate(t_high,n)),'ko');
            mRate_n = [mRate_n 60*nanmean(spkRate(t,n))];
            mRate_h = [mRate_h 60*nanmean(spkRate(t_high,n))];
            mRate_l = [mRate_l 60*nanmean(spkRate(t_low,n))];
            
%             pause
            drawnow
        end
    end
    subplot(248)
    axis equal; axis([0 30 0 30])
    line(xlim, ylim)
%     subplot(244)
%     axis equal; axis([0 0.30 0 0.30])
%     line(xlim, ylim)
%     subplot(243)
%     axis equal; axis([0 2 0 2])
%     line(xlim, ylim)
end