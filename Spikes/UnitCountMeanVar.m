function [mm,vv,BinDur] = UnitCountMeanVar(u, BinDur )
% UnitCountMeanVar calculates mean and variance of spike counts for a unit
%
% [mm,vv,BinDur] = UnitCountMeanVar(u ) calculates mean mm and variance vv for
% unit u. It does this stimulus by stimulus (bin duration is equal to
% stimulus duration)
%
% [mm,vv] = UnitCountMeanVar(u, BinDur ) lets you specify the bin duration
% in seconds. DEFAULT: [], which means take the whole stimulus.
%
% (assumes all stimuli have similar duration)
% 
% 2012-02 Matteo Carandini

if nargin < 2
    BinDur = [];
end

StimDur = median( u.stimdurs(:) ); % assumes all stimuli have similar duration

if isempty(BinDur)
    BinDur = StimDur;
end

nBins = floor(StimDur/BinDur);
BinTimes = 0:BinDur:(nBins*BinDur);

StimCounts = zeros( nBins, u.nstim, u.nrepeats);
for istim = 1:u.nstim
    for irepeat = 1:u.nrepeats
        Counts = histc( u.spiketimes{istim,irepeat}, BinTimes );
        StimCounts(:,istim,irepeat) = Counts(1:end-1);    
    end
end

     
mm = mean(StimCounts,    3);
vv =  var(StimCounts,[], 3);