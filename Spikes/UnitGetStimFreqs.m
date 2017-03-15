function fs = UnitGetStimFreqs(protocol, unit, stimNum)
% UnitGetStimFreqs returns the frequency of the stimulus
%
% [fs] = UnitGetStimFreqs(protocol, unit, stimNum)
% inputs are a protocol struct as obtained from ProtocolLoad, unit structure and an
% index for the stimulus. stimNum is typically 1 (e.g., if you are
% presenting on stimulus at a time), but can also be 2 or higher (e.g., if
% you are interested in the frequency of the 2nd orientation component of a
% plaid stimulus)
%
% 2008-07 LB
% 2009-03 MC bug fix: fields 'ncycles' and 'pfilefreqs' may not be populated
% 2010-06 MC made it deal with empty data files
% 2011-02 MC made it do more reasonable things before less reasonable ones

fs = [];

if ~isempty(protocol.estfreqs)
    fs = protocol.estfreqs;
end

if isempty(fs) && ...
        isfield(protocol, 'pfilefreqs') && ...
        size(protocol.pfilefreqs,1)>=stimNum && ...
        ~isempty(protocol.pfilefreqs(stimNum,:))
    fs = protocol.pfilefreqs(stimNum,:);
end

if isempty(fs) && ...
        isfield(protocol, 'ncycles') && ...
        size(protocol.ncycles,1)>=stimNum && ...
        ~isempty(protocol.ncycles(stimNum,:)) && ...
        ~any(isnan(protocol.ncycles(stimNum,:)))
    
    durs = unit.stimdurs;
    durs(durs==0)=NaN; % deal with empty files
    durs_data = nanmean(durs,2);
    fs = protocol.ncycles(stimNum,:)'./durs_data(:);
end

if isempty(fs)
    % give up, there is nothing we can do
    disp('WARNING !! Wasn`t able to find the temporal frequencies for this experiment');
    disp('Arbitrarly assuming the number of cycles to be 1 for each stimulus');
    
    durs = unit.stimdurs;
    durs(durs==0)=NaN; % deal with empty files
    durs_data = nanmean(durs,2);
    
    fs = (durs_data(:)*0+1)./durs_data(:);
end

if length(fs) ~= protocol.nstim
    error('<UnitGetStimFreqs> The protocol does not specify the frequency of the stimuli. Perhaps it was not periodic?');
end