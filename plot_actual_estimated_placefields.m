plot_idx = 1;
thres = 0.1;
figure('Position', [100 100 1600 1000])
for idx = 1:8
    list = find((nanmean(Posterior_all(idx).decoder.model.EV)>thres) | ...
        (nanmean(Posterior_all(idx).decoder_norm_decPos.model.EV)>thres));
    for icell = list %1:size(Posterior_all(idx).decoder.model.bestModel,1)
        tuning_actPos = 60*Posterior_all(idx).decoder.model.bestModel(icell,:);
        tuning_estPos = 60*Posterior_all(idx).decoder_norm_decPos.model.bestModel(icell,:);
        EV_actPos = round(100*nanmean(Posterior_all(idx).decoder.model.EV(:,icell)));
        EV_estPos = round(100*nanmean(Posterior_all(idx).decoder_norm_decPos.model.EV(:,icell)));
        
        if plot_idx>25
            drawnow; pause(2)
            figure('Position', [100 100 1600 1000])
            plot_idx = 1;
        end
        subplot(5,5,plot_idx)
        plot(1:2:100, tuning_actPos, 'k', 1:2:100, tuning_estPos,'r');
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        legend(num2str(EV_actPos), num2str(EV_estPos)...
            );%,'Position', 'best');
        plot_idx = plot_idx + 1;
        display(num2str([idx icell plot_idx]));
    end
end

