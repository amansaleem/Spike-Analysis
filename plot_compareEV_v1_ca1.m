function plot_compareEV_v1_ca1(maps)

thres = 0.01;
plot_each_exp = 0;

for iexp = 1:length(maps)
    t_V = (maps(iexp).V.EV_P>thres) | ...
        (maps(iexp).V.EV_S>thres) | ...
        (maps(iexp).V.EV_PS>thres);
    t_C = (maps(iexp).C.EV_P>thres) | ...
        (maps(iexp).C.EV_S>thres) | ...
        (maps(iexp).C.EV_PS>thres);
    
    %V1 mean EV
    meanEVs_P(iexp)     = nanmean(maps(iexp).V.EV_P(t_V));
    meanEVs_P_sem(iexp) = nansem(maps(iexp).V.EV_P(t_V));
    meanEVs_P_sig(iexp) = signrank(maps(iexp).V.EV_P(t_V));
    
    meanEVs_S(iexp)     = nanmean(maps(iexp).V.EV_S(t_V));
    meanEVs_S_sem(iexp) = nansem(maps(iexp).V.EV_S(t_V));
    meanEVs_S_sig(iexp) = signrank(maps(iexp).V.EV_S(t_V));
    
    meanEVs_PS(iexp)     = nanmean(maps(iexp).V.EV_PS(t_V));
    meanEVs_PS_sem(iexp) = nansem(maps(iexp).V.EV_PS(t_V));
    meanEVs_PS_sig(iexp) = signrank(maps(iexp).V.EV_PS(t_V));
    
    %CA1 mean EV
    meanEVsC_P(iexp)     = nanmean(maps(iexp).C.EV_P(t_V));
    meanEVsC_P_sem(iexp) = nansem(maps(iexp).C.EV_P(t_V));
    if sum(~isnan(maps(iexp).C.EV_P(t_V)))>0
        meanEVsC_P_sig(iexp) = signrank(maps(iexp).C.EV_P(t_V));
    else
        meanEVsC_P_sig(iexp) = nan;
    end
    
    
    meanEVsC_S(iexp)     = nanmean(maps(iexp).C.EV_S(t_V));
    meanEVsC_S_sem(iexp) = nansem(maps(iexp).C.EV_S(t_V));
    if sum(~isnan(maps(iexp).C.EV_S(t_V)))
        meanEVsC_S_sig(iexp) = signrank(maps(iexp).C.EV_S(t_V));
    else
        meanEVsC_S_sig(iexp) = nan;
    end
    
    meanEVsC_PS(iexp)     = nanmean(maps(iexp).C.EV_PS(t_V));
    meanEVsC_PS_sem(iexp) = nansem(maps(iexp).C.EV_PS(t_V));
    if sum(~isnan(maps(iexp).C.EV_PS(t_V)))
        meanEVsC_PS_sig(iexp) = signrank(maps(iexp).C.EV_PS(t_V));
    else
        meanEVsC_PS_sig(iexp) = nan;
    end
    
    meanDiffs_P_S(iexp) = nanmean((maps(iexp).V.EV_P(t_V))-(maps(iexp).V.EV_S(t_V)));
    meanDiffs_P_S_sem(iexp) = nansem((maps(iexp).V.EV_P(t_V))-(maps(iexp).V.EV_S(t_V)));
    meanDiffs_P_S_sig(iexp) = signrank((maps(iexp).V.EV_P(t_V))-(maps(iexp).V.EV_S(t_V)));
    
    meanDiffs_PS_P(iexp) = nanmean((maps(iexp).V.EV_PS(t_V))-(maps(iexp).V.EV_P(t_V)));
    meanDiffs_PS_P(iexp) = nansem((maps(iexp).V.EV_PS(t_V))-(maps(iexp).V.EV_P(t_V)));
    meanDiffs_PS_P_sig(iexp) = signrank((maps(iexp).V.EV_PS(t_V))-(maps(iexp).V.EV_P(t_V)));
    
    meanDiffs_PS_S(iexp) = nanmean((maps(iexp).V.EV_PS(t_V))-(maps(iexp).V.EV_S(t_V)));
    meanDiffs_PS_S_sem(iexp) = nansem((maps(iexp).V.EV_PS(t_V))-(maps(iexp).V.EV_S(t_V)));
    meanDiffs_PS_S_sig(iexp) = signrank((maps(iexp).V.EV_PS(t_V))-(maps(iexp).V.EV_S(t_V)));
    
    %CA1
    meanDiffsC_P_S(iexp) = nanmean((maps(iexp).C.EV_P(t_C))-(maps(iexp).C.EV_S(t_C)));
    meanDiffsC_P_S_sem(iexp) = nansem((maps(iexp).C.EV_P(t_C))-(maps(iexp).C.EV_S(t_C)));
    meanDiffsC_P_S_sig(iexp) = signrank((maps(iexp).C.EV_P(t_C))-(maps(iexp).C.EV_S(t_C)));
    
    meanDiffsC_PS_P(iexp) = nanmean((maps(iexp).C.EV_PS(t_C))-(maps(iexp).C.EV_P(t_C)));
    meanDiffsC_PS_P_sem(iexp) = nansem((maps(iexp).C.EV_PS(t_C))-(maps(iexp).C.EV_P(t_C)));
    meanDiffsC_PS_P_sig(iexp) = signrank((maps(iexp).C.EV_P(t_C))-(maps(iexp).C.EV_PS(t_C)));
    
    meanDiffsC_PS_S(iexp) = nanmean((maps(iexp).C.EV_PS(t_C))-(maps(iexp).C.EV_S(t_C)));
    meanDiffsC_PS_S_sem(iexp) = nansem((maps(iexp).C.EV_PS(t_C))-(maps(iexp).C.EV_S(t_C)));
    meanDiffsC_PS_S_sig(iexp) = signrank((maps(iexp).C.EV_S(t_C))-(maps(iexp).C.EV_PS(t_C)));
    
    if plot_each_exp
        figure;
        subplot(231)
        plot(maps(iexp).V.EV_P(t_V), maps(iexp).V.EV_S(t_V), 'ko')
        xlabel('V:EV_P'); ylabel('V:EV_S');
        title(['Experiment: ' num2str(iexp)]);
        subplot(232)
        plot(maps(iexp).V.EV_P(t_V), maps(iexp).V.EV_PS(t_V), 'ko')
        xlabel('V:EV_P'); ylabel('V:EV_{PS}');
        subplot(233)
        plot(maps(iexp).V.EV_S(t_V), maps(iexp).V.EV_PS(t_V), 'ko')
        xlabel('V:EV_S'); ylabel('V:EV_{PS}');
        
        subplot(234)
        plot(maps(iexp).C.EV_P(t_C), maps(iexp).C.EV_S(t_C), 'ko')
        xlabel('EV_P'); ylabel('EV_S');
        subplot(235)
        plot(maps(iexp).C.EV_P(t_C), maps(iexp).C.EV_PS(t_C), 'ko')
        xlabel('EV_P'); ylabel('EV_{PS}');
        subplot(236)
        plot(maps(iexp).C.EV_S(t_C), maps(iexp).C.EV_PS(t_C), 'ko')
        xlabel('EV_S'); ylabel('EV_{PS}');
        
        for n = 1:6
            subplot(2,3,n)
            set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
            axis equal;
            temp = max(max(xlim),max(ylim));
            axis([0 temp 0 temp])
            %     line([0 1], [1 0], 'linestyle','--', 'color','k');
            line(xlim, ylim, 'linestyle','--', 'color','k');
        end
        
        figure;
        subplot(231)
        hist((maps(iexp).V.EV_P(t_V))-..../(maps(iexp).C.EV_PS(:,t_C)),...
            (maps(iexp).V.EV_S(t_V)),..../(maps(iexp).C.EV_PS(:,t_C)),...
            30)
        title(['Mean: ' num2str(nanmean((maps(iexp).V.EV_P(t_V))-(maps(iexp).V.EV_S(t_V))))])
        ylabel([num2str(signrank((maps(iexp).V.EV_P(t_V))-(maps(iexp).V.EV_S(t_V))))]);
        
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('V:EV_P - V:EV_{S}')
        
        subplot(232)
        hist((maps(iexp).V.EV_PS(t_V))-..../(maps(iexp).C.EV_PS(:,t_C)),...
            (maps(iexp).V.EV_P(t_V)),..../(maps(iexp).C.EV_PS(:,t_C)),...
            30)
        title(['Mean: ' num2str(nanmean((maps(iexp).V.EV_PS(t_V))-(maps(iexp).V.EV_P(t_V))))])
        ylabel([ num2str(signrank((maps(iexp).V.EV_PS(t_V))-(maps(iexp).V.EV_P(t_V))))]);
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('V:EV_{PS} - V:EV_{P}')
        
        subplot(233)
        hist((maps(iexp).V.EV_PS(t_V))-..../(maps(iexp).C.EV_PS(:,t_C)),...
            (maps(iexp).V.EV_S(t_V)),..../(maps(iexp).C.EV_PS(:,t_C)),...
            30)
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('V:EV_{PS} - V:EV_{S}')
        title(['Mean: ' num2str(nanmean((maps(iexp).V.EV_PS(t_V))-(maps(iexp).V.EV_S(t_V))))])
        ylabel([ num2str(signrank((maps(iexp).V.EV_PS(t_V))-(maps(iexp).V.EV_S(t_V))))]);
        
        
        subplot(234)
        hist((maps(iexp).C.EV_P(t_C))-..../(maps(iexp).C.EV_PS(:,t_C)),...
            (maps(iexp).C.EV_S(t_C)),..../(maps(iexp).C.EV_PS(:,t_C)),...
            30)
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('EV_P - EV_{S}')
        title(['Mean: ' num2str(nanmean((maps(iexp).C.EV_P(t_C))-(maps(iexp).C.EV_S(t_C))))])
        ylabel([ num2str(signrank((maps(iexp).C.EV_P(t_C))-(maps(iexp).C.EV_S(t_C))))]);
        
        subplot(235)
        hist((maps(iexp).C.EV_PS(t_C))-..../(maps(iexp).C.EV_PS(:,t_C)),...
            (maps(iexp).C.EV_P(t_C)),..../(maps(iexp).C.EV_PS(:,t_C)),...
            30)
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('EV_{PS} - EV_{P}')
        title([num2str(nanmean((maps(iexp).C.EV_PS(t_C))-(maps(iexp).C.EV_P(t_C))))])
        ylabel([ num2str(signrank((maps(iexp).C.EV_P(t_C))-(maps(iexp).C.EV_PS(t_C))))]);
        
        subplot(236)
        hist((maps(iexp).C.EV_PS(t_C))-..../(maps(iexp).C.EV_PS(:,t_C)),...
            (maps(iexp).C.EV_S(t_C)),..../(maps(iexp).C.EV_PS(:,t_C)),...
            30)
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        xlabel('EV_{PS} - EV_{S}')
        title(['Mean: ' num2str(nanmean((maps(iexp).C.EV_PS(t_C))-(maps(iexp).C.EV_S(t_C))))])
        ylabel([ num2str(signrank((maps(iexp).C.EV_S(t_C))-(maps(iexp).C.EV_PS(t_C))))]);
        
        drawnow;
    end
end
%%
figure;
subplot(231)
errorbarxy(meanEVs_P,meanEVs_S,meanEVs_P_sem,meanEVs_S_sem,{'ko', 'k', 'k'})
xlabel('V:EV_P'); ylabel('V:EV_S');
title(num2str(signrank(meanDiffs_P_S)));
subplot(232)
errorbarxy(meanEVs_P,meanEVs_PS,meanEVs_P_sem,meanEVs_PS_sem,{'ko', 'k', 'k'})
xlabel('V:EV_P'); ylabel('V:EV_{PS}');
title(num2str(signrank(meanDiffs_PS_P)));
subplot(233)
errorbarxy(meanEVs_S,meanEVs_PS,meanEVs_S_sem,meanEVs_PS_sem,{'ko', 'k', 'k'})
xlabel('V:EV_S'); ylabel('V:EV_{PS}');
title(num2str(signrank(meanDiffs_PS_S)));

subplot(234)
errorbarxy(meanEVsC_P,meanEVsC_S,meanEVsC_P_sem,meanEVsC_S_sem,{'ko', 'k', 'k'})
xlabel('EV_P'); ylabel('EV_S');
title(num2str(signrank(meanDiffsC_P_S)));
subplot(235)
errorbarxy(meanEVsC_P,meanEVsC_PS,meanEVsC_P_sem,meanEVsC_PS_sem,{'ko', 'k', 'k'})
xlabel('EV_P'); ylabel('EV_{PS}');
title(num2str(signrank(meanDiffsC_PS_P)));
subplot(236)
errorbarxy(meanEVsC_S,meanEVsC_PS,meanEVsC_S_sem,meanEVsC_PS_sem,{'ko', 'k', 'k'})
xlabel('EV_S'); ylabel('EV_{PS}');
title(num2str(signrank(meanDiffsC_PS_S)));

for n = 1:6
    subplot(2,3,n)
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    axis equal;
    temp = max(max(xlim),max(ylim));
    axis([0 temp 0 temp])
    %     line([0 1], [1 0], 'linestyle','--', 'color','k');
    line(xlim, ylim, 'linestyle','--', 'color','k');
end
        