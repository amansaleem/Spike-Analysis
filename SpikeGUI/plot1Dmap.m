function plot1Dmap(input, icell, errorbaron, colorin, plot_index)

if nargin<3
    errorbaron = 1;
end
if nargin<4
    colorin = [0 0 0];
end
if nargin <5
    plot_index = gcf;
end

if ~isempty(input.model)
    tuningCurve = mean(input.model.tuning(icell).respModel,1).*input.sampleRate;
    errorset    = std(input.model.tuning(icell).respModel).*input.sampleRate;
    
    if errorbaron
        errorbar(plot_index,input.bins,tuningCurve,errorset,'color',colorin);
    else
        plot(plot_index,input.bins,tuningCurve,'color',colorin);
    end
    axis(plot_index,'tight');
    lims = get(plot_index, 'YLim');
    set(plot_index, 'YLim', [0 lims(2)]);
end