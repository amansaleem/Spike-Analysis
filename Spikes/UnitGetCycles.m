function [R, periods, dt] = UnitGetCycles( unit, protocol, resolution, options)
% UnitGetCycles the cycle averages of responses of a unit
%
% R = UnitGetCycles( unit, protocol ) returns the cycles of the responses
% 
% R = UnitGetCycles( unit, protocol, resolution ) let's you specify
% the resolution in seconds (DEFAULT: 0.025). Can be a number (the same
% across stimuli) or a vector (one per each stimulus).
%
% UnitGetCycles( unit, protocol, resolution, OPTS)
% lets you set the options OPTS
% OPTS(1) = 1 (DEFAULT: 0) Returns individual repeats. Useful if you want
% to fit entire sample set instead of sample set mean.
% OPTS(2) = 1 (DEFAULT: 0) Forces all stimuli to have the same number of
% bins.
%
% [R, periods] = UnitGetCycles( ... ) returns the periods as well
%
% [R, periods, dt] = UnitGetCycles( ... )
% gives the duration of the samples for each cycle (in seconds)
%
% part of Spikes
%
% 2000-12 MC
% 2001-03 MC figures out frequencies, returns periods
% 2001-08 VM added screenspec as an input, to compute true temp freq
%	added options
% 2001-10-10 VM now uses protocol.estfreqs instead of screenspec
% 2001-04-20 VM added times as an output, corrected bug of last bin
%  and resolution can vary across stimuli
% 2003-03 VM made it fit for traces, instead of times gives dt as an output
% 2003-12 MC made small modifications in resolution to yield integer n of bins
% 2004-02 VM corrected small bug in size of resolution
% 2004-02 MC added option to force resolution to be the same for all stimuli
% 2008-06 LB added an error message if spiketimes have the wrong format
%   (which is the case with a version of units created from cerebus)
% 2010-08 MC cleaned up, made it call UnitGetStimFreqs

% unit = UnitLoad( DIRS.spikes, 'CATZ008', 5, 9, 1, 2);
% protocol = ...

% Decide what datatype to use
if isfield(unit,'spiketimes') && ~isempty(unit.spiketimes)
   datatype = 'spiketimes';
elseif isfield(unit,'traces') && ~isempty(unit.traces)
   datatype = 'traces';
end

if nargin < 4 || isempty(options)
   options = [0 0];
end

if length(options)<2
    options(2) = 0;
end

if nargin < 3 || isempty(resolution)
   resolution = 0.010;
end

if ~isstruct(protocol)
   error('Second argument to UnitGetCycles should be a structure "protocol". You may be using old code.');
end

if isempty(protocol)
    error('Cannot work with an empty protocol');
end

if isfield(unit, 'spiketimes')
    try 
        [unit.spiketimes{1,:}]; %#ok<VUNUS>
    catch %#ok<CTCH>
        error('<UnitGetCycles> Units are in the wrong format. Run UnitGetCerebus again to fix it.');
    end
end

%%

data = getfield(unit,datatype);
nstim = size(data,1);
nrpts	= size(data,2);

if ~all( size(unit.stimdurs)==[nstim, nrpts] )
   error('UnitGetRates cannot work because data and stimdurs have different dimensions');
end


%% figure out the stimulus frequency for each stim

StimFreq = UnitGetStimFreqs(protocol, unit, 1);

%%

periods = 1./StimFreq;

if any(resolution < 0)
   error('Cannot have resolution less than 0 seconds');
end

if length(resolution) == 1
   resolutions = ones(size(periods))*resolution;
else
    resolutions = resolution;
end

% adjust the resolutions slightly so there is an integer number of bins
if options(2)
    nbins = round(median(periods./resolutions))*ones(size(periods));
else
    nbins = round(periods./resolutions);
end
resolutions = periods./nbins;

switch datatype
case 'spiketimes'
   if options(1)		% consider individual repeats
      R = cell(nstim,nrpts);
      dt = zeros(nstim,nrpts);
      for istim = 1:nstim
         for irpt = 1:nrpts
            dur = unit.stimdurs(istim,irpt);
            allspiketimes = rem( [unit.spiketimes{istim,irpt}], periods(istim));
            ncycles = dur/periods(istim);
            if isempty(allspiketimes)
               Nspikes = zeros(1, nbins(istim)+1);
            else
               % histc is fast because it is a Matlab mex function
               Nspikes = histc( allspiketimes/resolutions(istim), 0:nbins(istim) );
            end
            
            R{istim,irpt} = Nspikes(1:end-1) / (resolutions(istim) * ncycles); 
            
            % The sample duration
            dt(istim,irpt) = resolutions(istim);
         end
      end
   else
      R = cell(nstim,1);
      dt = resolutions;
      for istim = 1:nstim
         maxdur = max(unit.stimdurs(istim,:));
         allspiketimes = rem( [unit.spiketimes{istim,:}], periods(istim));
         ncycles = maxdur/periods(istim);
         if isempty(allspiketimes)
            Nspikes = zeros(1, nbins(istim)+1);
         else
            % histc is fast because it is a Matlab mex function
            Nspikes = histc( allspiketimes/resolutions(istim), 0:nbins(istim) );
         end
         
         R{istim} = Nspikes(1:end-1) / (resolutions(istim)*nrpts*ncycles);
         % removing the last one because it is spurious
         %if mod(period/resolutions(istim),1) > 0
         %   R{istim}(end) = R{istim}(end) / mod(period/resolutions(istim),1);
         %end
         % Spikes occur only in the first part of the last bin
      end
   end
   
case 'traces'
   
   cycleresp = cell(nstim,nrpts);
   dt = zeros(nstim,nrpts);
   period = 1./StimFreq'; % this variable already exists, it is called "periods"...
   
   % Choose the duration of a sample to be a fraction of a cycle
   NCYCLESAMPLEMIN = 10;
   ncyclesamples = max(ceil(period./resolutions),NCYCLESAMPLEMIN);
   cydt = period./ncyclesamples;
   
   for istim = 1:nstim
      for irpt = 1:nrpts
         tracedur = length(unit.traces{istim,irpt})*unit.sampledur;
         ncycles = ceil(tracedur/period(istim));
         
         % Find the time axis
         tracetime = (0:length(unit.traces{istim,irpt})-1)*unit.sampledur;
         cycletime = (0:ncyclesamples(istim)*ncycles-1)*cydt(istim);
         
         % Interpolate the trace
         traceinterp = interp1(tracetime,unit.traces{istim,irpt},cycletime);
         
         % Average over cycles
         cycleindex = reshape(1:ncyclesamples(istim)*ncycles,[ncyclesamples(istim) ncycles])';
         cycleresp{istim,irpt} = nanmean(traceinterp(cycleindex));

         % Replace the NaNs with zeros
         cycleresp{istim,irpt}(isnan(cycleresp{istim,irpt})) = 0;
      end
   end
   
   % Decide whether to average over repeats
   if options(1)
      R = cycleresp;
      dt = meshgrid(cydt,1:nrpts)';
   else
      dt = cydt;
      % Average over repeats
      R = cell(nstim,1);
      for istim = 1:nstim
         if nrpts > 1
            Rrpt = zeros(ncyclesamples(istim),nrpts);
            for irpt = 1:nrpts
               Rrpt(:,irpt) = cycleresp{istim,irpt};
               % deal with empty data files
               if unit.stimdurs(istim,irpt)==0
                   Rrpt(:,irpt) = NaN;
               end
            end
            R{istim} = nansum(Rrpt)/nrpts;
         else
            R{istim} = cycleresp{istim,1};
         end
      end
   end
end

   
