function [tuning]           = guantileTuning(condition, variable, zua, discreet_steps)
% function to quantile bin the data and get the tuning

if nargin<4
    discreet_steps = 12;
end

% getting an even distribution of bins
[distr bin_starts bin_centres] = eHist(variable(condition),discreet_steps);

bin_starts(diff(bin_centres)==0) = [];
distr(diff(bin_centres)==0) = [];
bin_centres(diff(bin_centres)==0) = [];

% finding the firing rates in each of the bins
if length(bin_starts)>2
    for n = 2:length(bin_starts)
        if bin_starts(n-1)==bin_starts(n)
            FR_std(n-1) =  nanstd(zua(variable==bin_starts(n-1) & condition));
            FR_sem(n-1) =  nansem(zua(variable==bin_starts(n-1) & condition));
            FR(n-1) = nanmean(zua(variable==bin_starts(n-1) & condition));
        else
            FR_std(n-1) =  nanstd(zua(variable>bin_starts(n-1) & variable<=bin_starts(n) & condition));
            FR_sem(n-1) =  nansem(zua(variable>bin_starts(n-1) & variable<=bin_starts(n) & condition));
            FR(n-1) = nanmean(zua(variable>bin_starts(n-1) & variable<=bin_starts(n) & condition));
        end
    end
    nFR = 60.*(FR);
    nFR_sem = (60.*FR_sem);
    nFR_std = (60.*FR_std);
    
    tuning.nFR = nFR;
    tuning.nFR_sem = nFR_sem;
    tuning.bin_starts = bin_centres;
else
    tuning.nFR = NaN;
    tuning.nFR_sem = NaN;
    tuning.nFR_starts = NaN;
end
end