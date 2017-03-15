function [R, E, resolution, dt] = UnitGetRates( unit, resolution, windowtype, dt, semflag, options)
% UnitGetRates the firing rate of responses in a Unit data structure
%
% [R, E] = UnitGetRates( unit ) returns firing rate R averaged across
% repeats and standard error E.
%
% [ ... ] = UnitGetRates( unit, resolution )
% lets you specify the resolution in seconds (default = 0.025 s).
%
% [ ... ] = UnitGetRates( unit, [], windowtype )
% lets you specify the window type:
% 'gauss' uses a gaussian window rather than a box (slow)
% 'active' uses active filtering rather than a box (faster)
% 'box' uses the classic box window (DEFAULT)
%
% [ ... ] = UnitGetRates( unit, resolution, 'gauss', dt )
% useful only for windowtype 'gauss':
% resolution is the sigma of the gaussian, dt is the sampling interval
% (set it to [] for other kinds of windowtype)
% DEFAULT if used is 1 ms (0.001)
%
% [ ... ] = UnitGetRates( unit, [], [], [], semflag )
% error E given as 'sem' (default) or 'std'
%
% [ ... ] = UnitGetRates( unit, resolution, windowtype, dt, semflag, options )
% lets you choose these options: [ KeepRepeats, SmartLength ]
% KeepRepeats: 0 to average across repeats [DEFAULT], 1 to keep them
% SmartLength: 0 to do something dumb, 1 for smart histogram length [DEFAULT]
%
% [R, E, resolution, dt] = UnitGetRates( ... )
% tells you the resolution of the histograms (trivial for windowtype 'box'
% and 'active') and the sampling interval dt.
%
% IF YOU NEED COUNTS INSTEAD OF RATES do something like this:
%  y = UnitGetRates(unit,0.001,[],[],[],[1 0]);
% for istim = 1:u.nstim, y{istim} = y{istim}/1000; end
%
% part of Spikes

% 2000-12 MC created
% 2001-10 VM computes sem for windowtype = 'box'
% 2001-12 MC added windowtype 'active' and computes sem for all window types
% 2002-03 MC fixed a little bug in the business of the gaussian window
% 2002-04 VM added times for windowtype 'box'
% 2003-03 VM now can use also traces. Replaced times in the output with dt.
% 2003-03 VM made it fit for traces
% 2003-04 VM small changes for windowtype 'gauss' case
% 2003-09 VM added dt as an input (for windowtype 'gauss')
% 2003-09 VM added semflag and options
% 2004-02 VM corrected bug in semflag (would never give you std)
% 2004-03 MC cleaned up the help text and rearranged the first few lines of code.
% 2004-12 VB now last bin of histogram also dropped when options=1
% 2006-01 VM added second option
% 2009-09 LB added error checking (length of unit struct)
% 2009-11 MC + ND made smart filtering the default, streamlined flag naming
% 2010-06 MC made it deal with missing files 
% 2011-03 ND put in a hack to shift spikes that happen at time 0 forward by
%            one sample spiketime processing. was breaking ungracefully. 
%            see below around line 200.
% 2013-08 ND for gaussian windowing changed filter to filtfilt so that
%            doesn't introduce a phase delay. see lines 211-212
% 2013-09 ND changed size checking for nstims and nrepeats to be the same
%            as unit.spiketimes rather than unit.stimdurs

% unit = UnitLoad( DIRS.spikes, 'CATZ008', 5, 9, 1, 2);

if nargin<2 || isempty(resolution)
    resolution = 0.025;
end

if nargin<3 || isempty(windowtype)
    windowtype = 'box';
end

if nargin < 4 || isempty(dt)
    dt = 0.001;
end

if nargin < 5 || isempty(semflag)
    semflag = 'sem';
end

if nargin < 6 || isempty(options)
    options = [0 1];
end

if length(options) == 1
    options(2) = 1;
end

if resolution<0
    error('Cannot have resolution less than 0 seconds');
end

if length(unit) > 1
    error('Only one instance of unit accepted for input');
end

%%

KeepRepeats = options(1);
SmartLength = options(2);

%%

% Decide whether to use spiketimes or traces
if isfield(unit,'spiketimes') && ~isempty(unit.spiketimes)
    datatype = 'spiketimes';
elseif isfield(unit,'traces') && ~isempty(unit.traces)
    datatype = 'traces';
end

if strcmp(datatype,'traces')
    resolution = [];
end

% data = getfield(unit,datatype);
data = unit.(datatype); % 2008-02-21

nstim = size(data,1);
nrpts	= size(data,2);

if ~all( size(unit.spiketimes)==[nstim, nrpts] ) % changed by NTD on 2013-09-11 should make no difference in most cases
% if ~all( size(unit.stimdurs)==[nstim, nrpts] )
    error('UnitGetRates cannot work because spiketimes and stimdurs have different dimensions');
end

R = cell(nstim,1);
E = cell(nstim,1);
% Rrpt = cell(nstim,1);  % 2008-02-21

switch datatype
    
    case 'spiketimes'
        
        switch windowtype
            case 'box'
                % The binsize
                dt = resolution;
                
                for istim = 1:nstim
                    maxdur = max(unit.stimdurs(istim,:));
                    
                    Rrpt = zeros(nrpts, ceil(maxdur/resolution)+1);
                    for irpt = 1:nrpts
                        spiketimes = [unit.spiketimes{istim,irpt}];
                        if ~isempty(spiketimes)
                            Rrpt(irpt,:) = histc( spiketimes/resolution, 0:ceil(maxdur/resolution) ) / resolution;
                        end
                        % handle missing files
                        if unit.stimdurs(istim,irpt)==0, Rrpt(irpt,:) = NaN; end
                    end
                    
                    if KeepRepeats
                        R{istim} = Rrpt;
                        E{istim} = [];
                    else
                        if nrpts > 1
                            switch semflag
                                case 'sem'
                                    R{istim} = nanmean(Rrpt) ; % MC 2011-05 R{istim} = nansum(Rrpt) / nrpts;
                                    E{istim} = nansem(Rrpt);
                                case 'std'
                                    R{istim} = nanmean(Rrpt) ; % MC 2011-05 R{istim} = nansum(Rrpt) / nrpts;
                                    E{istim} = nanstd(Rrpt);
                            end
                        else
                            R{istim} = Rrpt;
                            E{istim} = Rrpt * 0;
                        end
                    end % if KeepRepeats
                    
                    % The last bin contains spikes that fall EXACTLY on the last bin edge.
                    % Discard it as will never happen.
                    R{istim} = R{istim}(:,1:end-1);
                    E{istim} = E{istim}(:,1:end-1);
                    
                end % istim
                
            case 'gauss'
                % The time axis of the lowpass filter
                sigma = resolution;
                nfilt = ceil((sigma*4)/dt);
                ttfilt = (-nfilt:nfilt)*dt;
                
                % Make the filter
                myfilter = exp(-(ttfilt.^2/(2*(sigma)^2)))/sqrt(2*pi)/sigma;
                myfilter = myfilter/sum(myfilter);
                
                for istim = 1:nstim
                    maxdur = max(unit.stimdurs(istim,:));
                    
                    % Pick a good size for the histogram
                    if SmartLength
                        nsamples = ceil(maxdur/dt);
                        r = zeros(nsamples,nrpts);
                        Rrpt = zeros(nrpts,nsamples);
                    else
                        nsamples = ceil(maxdur/dt);
                        r = zeros(nsamples,nrpts);
                        Rrpt = zeros(nrpts,nsamples+nfilt*2);
                    end
                    
                    % Initialize the histogram
                    
                    % Find the times of each bin
                    for irpt = 1:nrpts
                        allspiketimes = ceil([unit.spiketimes{istim,irpt}]/dt);
                        % breaks ungracefully if a spike is at time 0. hack
                        % by shifting those ahead by one time sample
                        allspiketimes(allspiketimes==0) = 1;
                        r(allspiketimes,irpt) = 1/dt;
                        rrlong = [zeros(nfilt*2,1); r(:,irpt); zeros(nfilt*2,1)];
%                         rrfilt = filter(myfilter,1,rrlong);
                        rrfilt = filtfilt(myfilter,1,rrlong);
                        if SmartLength
                            Rrpt(irpt,:) = rrfilt(nfilt*3+1:end-nfilt)';
                        else
                            Rrpt(irpt,:) = rrfilt(nfilt*2+1:end)';
                        end
                        % handle missing files
                        if unit.stimdurs(istim,irpt)==0, Rrpt(irpt,:) = NaN; end
                    end
                    
                    
                    if KeepRepeats
                        R{istim} = Rrpt;
                        E{istim} = [];
                    else
                        if nrpts > 1
                            switch semflag
                                case 'sem'
                                    R{istim} = nanmean(Rrpt) ; % MC 2011-05 R{istim} = nansum(Rrpt) / nrpts;
                                    E{istim} = nansem(Rrpt);
                                case 'std'
                                    R{istim} = nanmean(Rrpt) ; % MC 2011-05 R{istim} = nansum(Rrpt) / nrpts;
                                    E{istim} = nanstd(Rrpt);
                            end
                        else
                            R{istim} = Rrpt;
                            E{istim} = Rrpt * 0;
                        end
                    end
                end
                
            case 'active'
                
                fc = 1/resolution;
                fs = fc*20;
                
                % The binsize
                dt = resolution;
                
                % fs = 1000; % sample rate
                % fc = 25; % cutoff rate
                
                for istim = 1:nstim
                    maxdur = max(unit.stimdurs(istim,:));
                    nsamples = ceil(maxdur*fs);
                    r = zeros(nrpts,nsamples);
                    Rrpt = zeros(nrpts,ceil(fc/fs*nsamples));
                    for irpt = 1:nrpts
                        SpikeTimes = unit.spiketimes{istim,irpt};
                        SpikeTimes(SpikeTimes<=0    ) = []; % added by MC 2011-06-15
                        SpikeTimes(SpikeTimes>maxdur) = []; % added by MC 2011-06-15
                        
                        if isnan(SpikeTimes)  % added by AAAS 2012-03-13
                            r(irpt,:) = NaN;
                        else
                            r(irpt,ceil(SpikeTimes*fs)) = fs;
                        end
                        Rrpt(irpt,:) = decimate(r(irpt,:),fs/fc);
                        % handle missing files
                        if unit.stimdurs(istim,irpt)==0, Rrpt(irpt,:) = NaN; end
                    end
                    
                    if KeepRepeats
                        R{istim} = Rrpt;
                        E{istim} = [];
                    else
                        if nrpts > 1
                            switch semflag
                                case 'sem'
                                    R{istim} = nanmean(Rrpt) ; % MC 2011-05 R{istim} = nansum(Rrpt) / nrpts;
                                    E{istim} = nansem(Rrpt);
                                case 'std'
                                    R{istim} = nanmean(Rrpt) ; % MC 2011-05 R{istim} = nansum(Rrpt) / nrpts;
                                    E{istim} = nanstd(Rrpt);
                            end
                        else
                            R{istim} = Rrpt;
                            E{istim} = Rrpt * 0;
                        end
                        E{istim} = max(0,E{istim});
                        R{istim} = max(0,R{istim});
                    end
                end
        end
        
    case 'traces'
        % The binsize
        dt = unit.sampledur;
        
        for istim = 1:nstim
            % Find the longest trace for this stimulus
            tracelength = zeros(nrpts,1);
            for irpt = 1:nrpts
                tracelength(irpt) = length(unit.traces{istim,irpt});
            end
            maxlength = max(tracelength);
            
            Rrpt = zeros(nrpts,maxlength);
            for irpt = 1:nrpts
                Rrpt(irpt,1:tracelength(irpt)) = unit.traces{istim,irpt};
                % handle missing files
                if unit.stimdurs(istim,irpt)==0, Rrpt(irpt,:) = NaN; end
            end
            
            if nrpts > 1
                switch semflag
                    case 'sem'
                        R{istim} = nanmean(Rrpt) ; % MC 2011-05 R{istim} = nansum(Rrpt) / nrpts;
                        E{istim} = nansem(Rrpt);
                    case 'std'
                        R{istim} = nanmean(Rrpt) ; % MC 2011-05 R{istim} = nansum(Rrpt) / nrpts;
                        E{istim} = nanstd(Rrpt);
                end
            else
                R{istim} = Rrpt;
                E{istim} = Rrpt * 0;
            end
        end
        
    otherwise
        error('Datatype unknown');
end


