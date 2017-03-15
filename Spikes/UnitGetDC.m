function [DCmean, DCerr] = UnitGetDC( unit, semflag, options )
% UnitGetDC returns the mean firing rate response 
%
% [DCmean, DCstd] = UnitGetDC( unit ) 
% returns the mean and the standard deviation of the responses
%
% [DCmean, DCsem] = UnitGetDC( unit, 'sem' ) 
% returns the standard error of the mean instead of the standard deviation
%
% DC = UnitGetDC( unit, [], 1)
% does not compute mean across repeats. Returns one response per repeat.
%
% This function is part of the Spikes toolbox.

% 2000-11 MC
% 2001-07 VB added nomeanflag
% 2001-09 VM changed nomeanflag into an option
% 2002-07 MC clarified the help text
% 2003-03 VM made it fit for traces
% 2008-02 LB handle repeats with missing start/stop time for cerebus units: substituted mean by nanmean
% 2010-06 MC handle missing data files
% 2010-06 cosmetic changes to deal with traces

% Decide what datatype to use
datatype = unit.datatype;

if nargin<3
   options = 0;
end

if nargin<2
   semflag = 'std';
end

data    = unit.(datatype); %getfield(unit,datatype);
nstim   = unit.nstim;
nrpts	= unit.nrepeats;

% if ~all(size(unit.stimdurs) == size(data))
%    error('UnitGetDC cannot work because data and stimdurs have different dimensions');
% end

% Compute the DC
switch datatype
    case 'spiketimes'
        spikecounts = zeros(nstim,nrpts);
        for istim = 1:nstim
            for irpt = 1:nrpts
                spikecounts(istim,irpt) = length(unit.spiketimes{istim,irpt});
            end
        end
        rr = spikecounts./unit.stimdurs;
    case 'traces'
        rr = zeros(nstim,nrpts);
        for istim = 1:nstim
            for irpt = 1:nrpts
                rr(istim,irpt) = mean(unit.traces{istim,irpt});
            end
        end
end

rr(unit.stimdurs==0) = NaN;

if options(1)
    DCmean =  rr; % one per repeat
else
    DCmean =  nanmean(rr,2); % one per stimulus (LB changed mean to nanmean 2008-Feb-06)
end


if nrpts>1 && ~options(1)
   switch semflag
   case 'std'
      DCerr =  nanstd(rr,0,2); % one per stimulus
   case 'sem'
      DCerr = nansem(rr')';
   otherwise
      error('Do not understand inputs to UnitGetDC');
   end
else
   DCerr = zeros(size(DCmean));
end
