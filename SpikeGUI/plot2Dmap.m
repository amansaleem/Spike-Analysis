function plot2Dmap(input, icell, plot_index)

if nargin <3
    plot_index = gca;
end

if ~isempty(input.model)
    tuningCurve = squeeze(mean(input.model.tuning(icell).respModel,1)).*input.sampleRate;
    
    axes(plot_index)
    imagesc(input.bins, input.binsB, tuningCurve');
    try
    set(gca, 'CLim',...
        [max([0 min(tuningCurve(:)) (nanmean(tuningCurve(:)) - 2*nanstd(tuningCurve(:)))])...
         min([max(tuningCurve(:)) (nanmean(tuningCurve(:)) + 2*nanstd(tuningCurve(:)))])]);
    catch
    end
    axis xy
    
end