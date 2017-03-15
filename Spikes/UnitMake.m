function unit = UnitMake(expt,chan,iprot,neighborhood)
% UnitMake creates a Unit data structure from Expt and Chan data structures
%
% UnitMake(expt,chan)
%
% part of Spikes
%
% 09/00 MC

ichan = chan.iexptchan;

icell = chan.cellids(iprot); % the cell id associated with the prototype

unit = struct(...
   'animal',expt.animal,...
   'ichan',expt.channames(ichan),...
   'icell',icell,...
   'iseries',expt.iseries,...
   'iexp',expt.iexp,...
   'stimdurs',expt.stimdurs);

unit.timestamp = expt.timestamp;

unit.prototype = chan.prots(iprot,:);
unit.neighborhood = neighborhood;

unit.spiketimes = cell(expt.nstim,expt.nrepeats);

if isempty(chan.cc), return; end

if isempty(find(chan.cc'==iprot)), return; end

if isempty(chan.spikes), return; end

protspikes = chan.spikes(find(chan.cc'==iprot));

if isempty(protspikes), return; end

for istim = 1:expt.nstim
   stimspikes = protspikes(find( [protspikes.istim]==istim ));
   if ~isempty(stimspikes)
      for irpt = 1:expt.nrepeats
         rptspikes = stimspikes(find( [stimspikes.irpt] == irpt));
         if ~isempty(rptspikes)
            unit.spiketimes{istim,irpt} = [rptspikes.t];
         end
      end
   end
end
