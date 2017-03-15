% SpikesExamples    Examples of batch scripts to analyze data
%
% 2009-08 MC revised
%
% %% EXAMPLE 1: PLOT THE TUNING
% 
% % Here is an example of a batch script that shows you the tuning for one experiment.
% % to run it, you have to add to your path the spikes toolbox and the matteobox toolbox.
% 
% SetDefaultDirs
% 
% animal = 'CATZ018';
% iseries = 1;
% iexp = 8;
% 
% protocol = ProtocolLoad( animal, iseries, iexp );
% 
% graphinfo = ProtocolGetGraphInfo(protocol);
% 
% units = UnitLoad( DIRS.spikes, animal, iseries, iexp );
% 
% for iunit = 1:length(units)
%    UnitPlotTuning( units(iunit), protocol, graphinfo );
% end
% 
% %% EXAMPLE 2: THE RAW DATA 
% 
% % Example of calling the most basic functions in the Spikes toolbox: load
% % the raw data and look at the spikes. This involves two types of
% % "objects": expt and chan
% % 
% % In general one does not need to look at these data, one only needs the
% % spike times. In essence one needs the objects "unit" and "protocol" see
% % example 1).
% 
% %---- basic information about the experiment
% animal = 'CAT001';
% iseries = 3;
% iexp = 1;
% 
% % ----------- get the raw data ----------------
% expt = ExptLoad(animal, iseries, iexp);
% 
% %---- look at one response
% istim = 1;
% irpt = 1;
% figure; plot(expt.data{istim}{irpt});
% 
% % ---------- get the analyzed data -------
% chan = ChanLoad(DIRS.spikes,animal,iseries,iexp);
% 
% % plot a prototype spike
% figure; plot(chan.tt,chan.prots(1,:));
