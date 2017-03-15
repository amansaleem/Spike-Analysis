function newchans = ChanImport(expt,chans)
% ChanImport import a new set of chans for an existing expt
%
%	newchans = ChanImport(expt,chans)
%  tries to import a new chans for an existing expt
%
% MC 00

for ichan = 1:expt.nchans
   channame = expt.channames(ichan);
   % see if there is a chan that matches this channame
   if isempty([chans.channame])
      % it is one of the old chans, there was no channame saved
      if length(chans)>=ichan
         newchans(ichan) = chans(ichan);
         newchans(ichan).channame = channame;
      else
         % this ichan does not match. may be new.
         fprintf(1,'Spike parameters are not present for channel %d\n', expt.channames(ichan) );  
         newchans(ichan) = ChanInitialize(expt,ichan);
      end
   else
      theichan = find([chans.channame] == channame);
      if length(theichan)==1
         oldchan = chans(theichan);
         newchans(ichan) = ChanInitialize;
         thefields = fieldnames(newchans(ichan));
         nfields = length(thefields);
         for ifield = 1:nfields
            thefield = thefields{ifield};
            newchans(ichan)= setfield(newchans(ichan),thefield,getfield(oldchan,thefield));
         end
         % newchans(ichan) = chans(theichan);
         newchans(ichan).iexptchan = ichan;
      else
         % this ichan does not match. may be new.
         fprintf(1,'Spike parameters are not present for channel %d\n', expt.channames(ichan) );  
         newchans(ichan) = ChanInitialize(expt,ichan);
      end
   end
end

