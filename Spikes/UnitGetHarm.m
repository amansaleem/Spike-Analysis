function [Hmean, Herr, fs] = UnitGetHarm( f, unit, protocol, semflag, broadflag, options, stimNum, fs)
% UnitGetHarm returns a frequency component of the firing rate response
%
% [Hmean, Hstd, fs] = UnitGetHarm( f, unit )
% returns the mean and the standard deviation of f-th harmonic
%
% [Hmean, Hstd, fs] = UnitGetHarm( f, unit, protocol ) lets you specify the
% protocol if you really want to (unnecessary). DEFAULT: [].
%
% For example, UnitGetHarm(1,...) returns the first harmonic, and UnitGetHarm(2,...) returns the second harmonic.
% Harmonics are reported PEAK to PEAK
%
% [Hmean, Hstd, fs] = UnitGetHarm( f, unit, protocol, 'std' )
% does the same
%
% [Hmean, Hsem, fs] = UnitGetHarm( f, unit, protocol, 'sem' )
% returns the mean and the standard error of the mean of f-th harmonic the responses
%
% [Hmean, Hsem, fs] = UnitGetHarm( f, unit, protocol, ..., 'broad' )
% looks around for the biggest response near the indicated frequencies. Useful for those
% case when you are not entirely sure of the stimulus frequency.
%
% [Hmean, Herr, fs] = UnitGetHarm( f, unit, protocol,[],[], 1)
% does not compute mean across repeats. returns one coefficient per repeat
% useful if you want to fit entire sample set instead of sample set mean
%
% [Hmean, Herr, fs] = UnitGetHarm( f, unit, protocol,'std','narrow', 0, stimNum)
% takes the frequency if the stimulus specified in stimNum - helpful if you
% are using more than one frequency (make sure that protocol.ncycles and
% protocol.pfilefreqs return an entry for each stimulus)
%
% [Hmean, Herr, fs] = UnitGetHarm( f, unit, protocol,'std','narrow', 0, [], fs)
% lets you specify the frequency of interest whose harmonics you are
% interested in - useful if you want to measure the amount of noise off from the
% freq of interest. fs must have the same dimensions as unit.spiketimes or
% unit.traces.
%
% part of Spikes
%
% 2000-11 Matteo Carandini (UnitGetF1)
% 2001-07 MC generalized to any harmonic
% 2000-07 VB broadflag can now take value 'nomean'
% 2000-08 VM added screenspec as an input and fs as an output
%            added options: options == 1 corresponds to broadflag == 'nomean' in the
%            previous version
% 2001-10-10 VM now uses protocol.estfreqs instead of screenspec
% 2002-03-04 VM & MC, integrates over an integer number of cycles, not over the
%                     entire stimulus duration
% 2003-03 VM made it fit for traces
% 2003-04 MC extended it so f can be zero (the mean)
% 2008-06 LB adjusted code for the use of 2 stim freqs and corrected some
% 	         Matlab warnings
% 2008-07 LB replaced code to determine stim freqs by the function
%            UnitGetStimFreqs
% 2008-07 LB added the input option fs
% 2011-02 MC added case with nargin = 2
% 2013-09 ND changed the line where we ignore stimuli of zero duration. now
%            limited to the known run repeats and stimuli
%
% unit = UnitLoad( DIRS.spikes, 'CATZ008', 5, 2, 1, 2);
% protocol = ProtocolLoad( unit.animal, unit.iseries, unit.iexp );

% Decide what datatype to use
datatype = unit.datatype;

if nargin < 8
    fs = [];
end
if nargin < 7
    stimNum = 1;
end

if nargin < 6
    options = 0;
end

if nargin < 5
    broadflag = 'narrow';
end

if nargin < 4
    semflag = 'std';
end

if nargin < 3
    protocol = [];
end

%%

if isempty(protocol)
    protocol = ProtocolLoad(unit);
end

%%
    
    
if f<0, error('f must be positive'); end

if f==0
    [Hmean, Herr] = UnitGetDC( unit, semflag, options );
    fs = [];
    return
end

%-------- figure out the stimulus frequency for each stim
if isempty(fs)
    fs = UnitGetStimFreqs(protocol,unit,stimNum);
end
%-------- update fs to be the frequency the user asked for
fs = f*fs;
%-------------------------

%data = getfield(unit,datatype);
nstim = size(unit.(datatype),1);
nrpts	= size(unit.(datatype),2);

if size(unit.stimdurs)~=size(unit.(datatype))
    error('data and stimdurs have different dimensions');
end

if strcmp(broadflag,'broad')
    if strcmp(datatype,'spiketimes')
        p = 0.01; % max error in estimation of actual period, in s
        R = UnitGetRates( unit, 0.005, 'gauss');
        newfs = fs; % for allocation and default value
        % loop over nonblank stimuli
        for istim = setdiff(1:nstim,protocol.blankstims)
            nsamples = length(R{istim});
            duration = nsamples/1000;
            absrft = abs(fft(R{istim})) /(nsamples/2);
            ff = freq(nsamples,duration);
            freqmin = 1/(1/fs(istim)+p);
            freqmax = 1/(1/fs(istim)-p);
            ii = find(ff>=freqmin & ff<=freqmax);
            if ~isempty(ii)
                % OLD: indmax = findmax(absrft(ii));
                [ mymax, indmax ] = max(absrft(ii));
                indmax = indmax(1);
                newfs(istim) = ff(ii(indmax));
            end
        end
        fs = newfs;
    else
        warning('UnitGetHarm:CannotDoBroad','Don`t know how to compute tuning in `broad` modus if datatype is not spiketimes');
    end
end

%%

switch datatype
    
    case 'spiketimes'
        re = zeros(nstim,nrpts);
        im = zeros(nstim,nrpts);
        for istim = 1:nstim
            for irpt = 1:nrpts
                tt = unit.spiketimes{istim,irpt};
                if f == 0
                    dur = unit.stimdurs(istim,irpt);
                else
                    period = 1/fs(istim); % in seconds
                    nperiods = floor(unit.stimdurs(istim,irpt)/period);
                    if nperiods < 1
                        dur = period;
                    else
                        dur = nperiods*period;
                    end
                    tt(tt>dur) = [];
                end
                re(istim,irpt) = sum(cos(2*pi*tt*fs(istim)))/dur;
                im(istim,irpt) = sum(sin(2*pi*tt*fs(istim)))/dur;
            end
        end
        H = 2 * (re+sqrt(-1)*im); % peak to peak
        
    case 'traces'
        % HACK: histogetharm -- called only once, is defined below.
        dt = unit.sampledur;
        H = nan(nstim,nrpts);
        for istim = 1:nstim
            for irpt = 1:nrpts
                tt = (0:length(unit.traces{istim,irpt})-1)*dt;
                if ~isempty(unit.traces{istim,irpt})
                    H(istim,irpt) = histogetharm(1,unit.traces{istim,irpt},tt,fs(istim));
                end
            end
        end
end

% H(unit.stimdurs == 0) = NaN;
H(unit.stimdurs(1:nstim,1:nrpts) == 0) = NaN; % changed by NTD
if options(1) % does not compute the mean over repeats
    Hmean = H;
else
    Hmean =  nanmean(H,2); % one per stimulus
end

if nrpts>1 && ~options(1)
    switch semflag
        case 'std'
            Herr =  nanstd(abs(H),0,2); % one per stimulus
        case 'sem'
            Herr = nansem(abs(H)')';
        otherwise
            error('Do not understand inputs to UnitGetHarm');
    end
else
    Herr = zeros(size(Hmean));
end

function h = histogetharm(harm,histo,tt,fs)
% histogetharm compute the power at the harmonic harm in the signal histo
% 
% THIS SHOULD NOT BE ITS OWN FUNCTION -- IT IS CALLED ONCE
% WAS PART OF VALERIO'S code
% 
% Inputs are:
% harm: the harmonics to compute
% histo: the signal
% tt: the time axis
% fs: the fundamental frequency 
%
% Assumes that all time bins have the sime length
%
% The output is a complex number (amplitude and phase)
% 
% 20 Feb 2003 Valerio Mante
% 
% h = histogetharm(harm,histo,tt,fs)

% The duration of a bin
dt = tt(2) - tt(1);

% The number of samples
nsamples = length(tt);

% The duration of the signal
dur = nsamples * dt;

% The frequency at which to compute the power
omega = harm*fs;

% Decide whether you can use to entire stimulus
if harm == 0
    % Use everything to compute the dc
    gooddur = dur;
    touse = 1:nsamples;
else
    % Use an entire number of periods to compute a harmonic
    period = 1/fs;
    nperiods = floor(dur/period);
    if nperiods < 1
        gooddur = period;
    else
        gooddur = nperiods*period;
    end
    touse = 1:round(gooddur/dt);
end
re = sum(histo(touse).*cos(2*pi*tt(touse)*omega)) * dt / gooddur;
im = sum(histo(touse).*sin(2*pi*tt(touse)*omega)) * dt / gooddur;
if harm > 0
    h = 2 * (re+sqrt(-1)*im);
else
    h = (re+sqrt(-1)*im); 
end



