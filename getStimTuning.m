function [tuning, stim, p, cellIDs] = getStimTuning(animal, iseries, iexp, shank_id, flag_plot, addInfo)

if nargin<5
    flag_plot = 1;
end
if nargin<4
    shank_id = 8;
end
if nargin<6
    addInfo = '1';
end

[outputMatrix, p, cellIDs] = getStimSpiketimes(animal, iseries, iexp, shank_id, addInfo);

stim = p.pars(p.activepars{1},:);

tuning.mean = zeros(length(outputMatrix),size(p.seqnums,1));
tuning.sem  = zeros(length(outputMatrix),size(p.seqnums,1));
tuning.active = zeros(1,length(outputMatrix));

for icell = 1:length(outputMatrix)
    if ~isempty(outputMatrix{icell})
        tuning.mean(icell,:) = mean(outputMatrix{icell}');
        tuning.sem(icell,:)  = sem(outputMatrix{icell}');
        tuning.active(icell) = 1;
    end
end

if flag_plot
    n = 1;
    for icell = 1:length(outputMatrix)
        if tuning.active(icell)==1 ...
            ;% && length(intersect(cellIDs(icell), [32 41 45 58 90 93 108 49 115 61 99 101 114 31]))>0 
            subplot(6,6,n)
            errorbar(stim(1:end-1), tuning.mean(icell,1:end-1),tuning.sem(icell,1:end-1),'bo-');
            hold on;
%             pars = fitori(stim(1:end-1), tuning.mean(icell,1:end-1));
%             plot(0:360, oritune(pars, 0:360), 'r-' );
            axis tight;
            line(xlim, [tuning.mean(icell,end) tuning.mean(icell,end)], 'color', 'k', 'linestyle','--');
            axis([xlim 0 max(ylim)]);
            set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%             set(gca, 'XTick', [0 120 240 360])
            title(num2str(cellIDs(icell)));
            hold off
            n = n + 1;
        end
    end
end