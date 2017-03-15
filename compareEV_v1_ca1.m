function maps = compareEV_v1_ca1(type, speedThres)

if nargin<1
    type = 'med';
end

if nargin<2
    speedThres = 5;
end
smthWin = 250;
box_filt = 0;
expt_list = getExptList;

for iexp = 1:length(expt_list)
    % Load files
    es = VRLoadMultipleExpts(expt_list(iexp).animal, expt_list(iexp).iseries,...
        expt_list(iexp).expt_list, 'SPIKES', [20 40], expt_list(iexp).shank_list_V1,'1');
    esC = VRLoadMultipleExpts(expt_list(iexp).animal, expt_list(iexp).iseries,...
        expt_list(iexp).expt_list, 'SPIKES', [20 40], expt_list(iexp).shank_list_CA1,'1');
    % For expts 1 and 2 get subset of cells
    if (iexp==1) | (iexp==2)
        es.spikeTrain = es.spikeTrain(:,expt_list(iexp).V1_cell_list);
        esC.spikeTrain = esC.spikeTrain(:,expt_list(iexp).CA1_cell_list);
    end
    % Remove interneurons
    esC.spikeTrain(:,expt_list(iexp).interneurons) = [];
    
    % Smooth spike trains
    for icell = 1:size(es.spikeTrain,2)
        es.spikeTrain(:,icell) = smthInTime(es.spikeTrain(:,icell), 60, smthWin);
    end
    for icell = 1:size(esC.spikeTrain,2)
        esC.spikeTrain(:,icell) = smthInTime(esC.spikeTrain(:,icell), 60, smthWin);
    end
    
    %     Get the conditions
    base = es.traj>1 & es.contrast~=0 & ~isnan(es.traj) & es.gain==1 & es.roomLength==1 ...
        & es.smthBallSpd>speedThres & es.trajspeed>=0;
    t_med      = round(es.contrast*100)/100==0.6 & base;
    t_low  = round(es.contrast*100)/100==0.18 & base;
    t_high = round(es.contrast*100)/100==0.72 & base;
    t_gray = es.traj~=0 & es.contrast==0 & ~isnan(es.traj) & es.smthBallSpd>5;
    
    switch type
        case 'low'
            t = t_low;
        case 'med'
            t = t_med;
        case 'high'
            t = t_high;
    end
    %     get the variables, position and speed
    X_P = es.traj;
    X_S = es.smthBallSpd;
    
    % Initialise the models
    mapP_V = oneDimMap;
    mapS_V = oneDimMap;
    mapPS_V = twoDimMap;
    
    mapP_C = oneDimMap;
    mapS_C = oneDimMap;
    mapPS_C = twoDimMap;
    
    Y  = es.spikeTrain;
    YC = esC.spikeTrain;
    
    % train the maps
    mapP_V = mapP_V.trainSpikeMap([X_P(t)], Y(t,:), 0);
    mapS_V = mapS_V.trainSpikeMap([X_S(t)], Y(t,:), 0);
    mapPS_V = mapPS_V.trainSpikeMap([X_P(t) X_S(t)], Y(t,:), 0);
    
    mapP_C = mapP_C.trainSpikeMap([X_P(t)], YC(t,:), 0);
    mapS_C = mapS_C.trainSpikeMap([X_S(t)], YC(t,:), 0);
    mapPS_C = mapPS_C.trainSpikeMap([X_P(t) X_S(t)], YC(t,:), 0);
    
    t_C = nanmean(mapPS_C.model.EV)>0.001;
    t_V = nanmean(mapPS_V.model.EV)>0.001;
    
    maps(iexp).t    = t;
    maps(iexp).info = expt_list(iexp);
    
    maps(iexp).V.P  = mapP_V;
    maps(iexp).V.S  = mapS_V;
    maps(iexp).V.PS = mapPS_V;
    maps(iexp).V.EV_PS = nanmean(mapPS_V.model.EV);
    maps(iexp).V.EV_P  = nanmean(mapP_V.model.EV);
    maps(iexp).V.EV_S  = nanmean(mapS_V.model.EV);
    
    maps(iexp).C.P  = mapP_C;
    maps(iexp).C.S  = mapS_C;
    maps(iexp).C.PS = mapPS_C;
    maps(iexp).C.EV = nanmean(mapPS_C.model.EV);
    maps(iexp).C.EV_PS = nanmean(mapPS_C.model.EV);
    maps(iexp).C.EV_P  = nanmean(mapP_C.model.EV);
    maps(iexp).C.EV_S  = nanmean(mapS_C.model.EV);
    
%     figure;
%     subplot(231)
%     plot(maps(iexp).V.EV_P, maps(iexp).V.EV_S, 'ko')
%     xlabel('V:EV_P'); ylabel('V:EV_S'); 
%     title(['Experiment: ' num2str(iexp)]);
%     subplot(232)
%     plot(maps(iexp).V.EV_P, maps(iexp).V.EV_PS, 'ko')
%     xlabel('V:EV_P'); ylabel('V:EV_{PS}'); 
%     subplot(233)
%     plot(maps(iexp).V.EV_S, maps(iexp).V.EV_PS, 'ko')
%     xlabel('V:EV_S'); ylabel('V:EV_{PS}'); 
%     
%     subplot(234)
%     plot(maps(iexp).C.EV_P, maps(iexp).C.EV_S, 'ko')
%     xlabel('EV_P'); ylabel('EV_S'); 
%     subplot(235)
%     plot(maps(iexp).C.EV_P, maps(iexp).C.EV_PS, 'ko')
%     xlabel('EV_P'); ylabel('EV_{PS}'); 
%     subplot(236)
%     plot(maps(iexp).C.EV_S, maps(iexp).C.EV_PS, 'ko')
%     xlabel('EV_S'); ylabel('EV_{PS}'); 
%     
%     for n = 1:6
%         subplot(2,3,n)
%         set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%         axis equal;
%         temp = max(max(xlim),max(ylim));
%         axis([0 temp 0 temp])
%         %     line([0 1], [1 0], 'linestyle','--', 'color','k');
%         line(xlim, ylim, 'linestyle','--', 'color','k');
%     end    
%     
%     figure;
%     subplot(231)
%     hist(nanmean(mapP_V.model.EV)-..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         nanmean(mapS_V.model.EV),..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         30)
%     title(['Mean: ' num2str(nanmean(nanmean(mapP_V.model.EV)-nanmean(mapS_V.model.EV)))])
%         ylabel([num2str(signrank(nanmean(mapP_V.model.EV)-nanmean(mapS_V.model.EV)))]);
%     
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%     xlabel('V:EV_P - V:EV_{S}')
%     subplot(232)
%     hist(nanmean(mapP_V.model.EV)-..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         nanmean(mapPS_V.model.EV),..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         30)
%     title(['Mean: ' num2str(nanmean(nanmean(mapP_V.model.EV)-nanmean(mapPS_V.model.EV)))])
%         ylabel([ num2str(signrank(nanmean(mapP_V.model.EV)-nanmean(mapPS_V.model.EV)))]);
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%     xlabel('V:EV_P - V:EV_{PS}')
%     subplot(233)
%     hist(nanmean(mapS_V.model.EV)-..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         nanmean(mapPS_V.model.EV),..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         30)
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%     xlabel('V:EV_S - V:EV_{PS}')
%     title(['Mean: ' num2str(nanmean(nanmean(mapS_V.model.EV)-nanmean(mapPS_V.model.EV)))])
%         ylabel([ num2str(signrank(nanmean(mapS_V.model.EV)-nanmean(mapPS_V.model.EV)))]);
%     
%     
%     subplot(234)
%     hist(nanmean(mapP_C.model.EV)-..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         nanmean(mapS_C.model.EV),..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         30)
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%     xlabel('EV_P - EV_{S}')
%     title(['Mean: ' num2str(nanmean(nanmean(mapP_C.model.EV)-nanmean(mapS_C.model.EV)))])
%         ylabel([ num2str(signrank(nanmean(mapP_C.model.EV)-nanmean(mapS_C.model.EV)))]);
%     
%     subplot(235)
%     hist(nanmean(mapP_C.model.EV)-..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         nanmean(mapPS_C.model.EV),..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         30)
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%     xlabel('EV_P - EV_{PS}')
%     title([num2str(nanmean(nanmean(mapP_C.model.EV)-nanmean(mapPS_C.model.EV)))])
%         ylabel([ num2str(signrank(nanmean(mapP_C.model.EV)-nanmean(mapPS_C.model.EV)))]);
%     
%     subplot(236)
%     hist(nanmean(mapS_C.model.EV)-..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         nanmean(mapPS_C.model.EV),..../nanmean(mapPS_C.model.EV(:,t_C)),...
%         30)
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%     xlabel('EV_S - EV_{PS}')
%     title(['Mean: ' num2str(nanmean(nanmean(mapS_C.model.EV)-nanmean(mapPS_C.model.EV)))])
%         ylabel([ num2str(signrank(nanmean(mapS_C.model.EV)-nanmean(mapPS_C.model.EV)))]);
%     ylabel(num2str(iexp));
%     
%     drawnow;
%     pause
    
end
    function expt_list = getExptList
        idx = 0;
        idx = idx + 1;
        expt_list(idx).animal     = 'M130920_BALL';
        expt_list(idx).iseries    = 1025;
        expt_list(idx).expt_list  = 102:103;
        expt_list(idx).shank_list_V1  = 0:7;
        expt_list(idx).shank_list_CA1 = 0:7;
        expt_list(idx).V1_cell_list = 153:249;
        expt_list(idx).CA1_cell_list = 1:152;
        expt_list(idx).interneurons = [17, 18,  9,11,12, 44, 45, 41, 78, 102, 104, 129, 113, 140, 145];
        
        idx = idx + 1;
        expt_list(idx).animal     = 'M130918_BALL';
        expt_list(idx).iseries    = 1030;
        expt_list(idx).expt_list  = 103:105;
        expt_list(idx).shank_list_V1  = 0:7;
        expt_list(idx).shank_list_CA1 = 0:7;
        expt_list(idx).V1_cell_list = 190:321;
        expt_list(idx).CA1_cell_list = 1:189;
        expt_list(idx).interneurons = [15, 20, 30, 39, 71, 132];
        
        idx = idx + 1;
        expt_list(idx).animal     = 'M140501_BALL';
        expt_list(idx).iseries    = 530;
        expt_list(idx).expt_list  = 104:106;
        expt_list(idx).shank_list_V1  = 8;
        expt_list(idx).shank_list_CA1 = 0:7;
        expt_list(idx).interneurons = [4 8 36 39 56 62 72 78 84 95 98 101 99 102 103 105 131 155 163 167 207 217 222 229 231 266 270];
        
        idx = idx + 1;
        expt_list(idx).animal     = 'M140501_BALL';
        expt_list(idx).iseries    = 531;
        expt_list(idx).expt_list  = 103:106;
        expt_list(idx).shank_list_V1  = 8;
        expt_list(idx).shank_list_CA1 = 0:7;
        expt_list(idx).interneurons = [11 15 16 44 53 86 94 98];
        
        idx = idx + 1;
        expt_list(idx).animal     = 'M140501_BALL';
        expt_list(idx).iseries    = 601;
        expt_list(idx).expt_list  = 103:106;
        expt_list(idx).shank_list_V1  = 8;
        expt_list(idx).shank_list_CA1 = 0:7;
        expt_list(idx).interneurons = [4 9 10 11 12 20 24 29 47 50 52 60 64 71 84];
        
        idx = idx + 1;
        expt_list(idx).animal     = 'M140501_BALL';
        expt_list(idx).iseries    = 602;
        expt_list(idx).expt_list  = 102:106;
        expt_list(idx).shank_list_V1  = 8;
        expt_list(idx).shank_list_CA1 = 0:7;
        expt_list(idx).interneurons = [5 8 10 18 24 29 35 38 40 52 54 59 60 61 77 90 92 93 107 109 127];
        
        idx = idx + 1;
        expt_list(idx).animal     = 'M140502_BALL';
        expt_list(idx).iseries    = 603;
        expt_list(idx).expt_list  = [107:108 110];
        expt_list(idx).shank_list_V1  = 8;
        expt_list(idx).shank_list_CA1 = 0:7;
        expt_list(idx).interneurons = [3 8 9 15 57 63 64 66 74 116 131 138 171 177 182 185 188];
        
        idx = idx + 1;
        expt_list(idx).animal     = 'M140502_BALL';
        expt_list(idx).iseries    = 604;
        expt_list(idx).expt_list  = 107:110;
        expt_list(idx).shank_list_V1  = 8;
        expt_list(idx).shank_list_CA1 = 0:7;
        expt_list(idx).interneurons = [8 22 25 26 19 31 40 52 53 66 91 97 109];
        
    end

end