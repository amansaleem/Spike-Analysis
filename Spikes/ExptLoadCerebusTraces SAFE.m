function expt = ExptLoadCerebusTraces(varargin)
% EXPTLOADCEREBUSTRACES: takes an ExptTag and returns an "expt" structure.
% This function is meant for continuous data (i.e. ns* files, not nev files!)
%
% Usage:
% Passes an ExptTag (with fields animal, iseries, iexp) as the argument:
% expt = ExptLoadCerebusTraces(ExptTag)
%
% Can also have as arguments the animal, iseries, iexp values themselves:
% expt = ExptLoadCerebusTraces(animal,iseries,iexp)
% 
% If you want to internally reference:
% expt = ExptLoadCerebusTraces(animal,iseries,iexp,refChan)
% NOTE that if you do this then you must use separate arguments
% for the animal, iseries and iexp and then a fourth argument is the
% channel to use as a reference. Note, this channel is the one after the
% probe has been unmapped, i.e., 3 means to reference from channel 3 - the
% third from the top in an A1x16 michigan probe.

%
% NOTICE THAT TRACES ARE DOWNSAMPLED BY A FACTOR OF 2 UPON LOADING

% MS 2010-03-22 created (adapted from UtahGetLFP)
% MS 2010-03-24 made some minor changes to calculating sample rate and stimulus durations
% MC 2010-06-14 changed the order of opening ns5 and nev (trouble with memmap)
% MC 2010-06-14 introduced a downsample by a factor of 2
% MC 2010-06-16 attempted to fix bug in timestamp (search for BUG)
% MC 2010-06-16 saves the data on zserver as a matlab file (an "expt")
% MC 2010-06-16 From then on, when you call it, it loads the mat file directly. 
% AS 2010-07-19 Added the functionality to load a Michigan map file
% ND 2010-09-01 Suggest we take out the call to ProtocolAbsentMatfiles
%               because we no longer collect data using Multispike and this leads to
%               warnings that really don't matter.
% ND 2010-10-10 took out call to ProtocolAbsentMatfiles
% ND 2010-10-21 allowed for Cerebus files to be named more appropriately
%               and for them to be saved by repeat rather than by experiment.
% MS 2010-10-25 added another input argument: StimTime, that is a vector of length 2 
%               (e.g. [3 1]; 3 s of activity before stimulus onset, and 1 s of activity after stimulus offset
% ND 2010-12-09 Fixed bugs introduced by AS 2010-07-19 (last site was not
%               remapped) and by ND 2010-10-21 (weird assignment of repeats).
% MS 2010-12-14 when StimTime exists (see 2010-10-25), this code is always run, even when an expt file already exists 
%               (however the new expt struct is not saved)
% AA 2011-02-01 deleted photodiode time stamps (PDtimestamps) that are closer than 5 ms
% MS 2011-02-23 added another input argument: refChan, that subtracts the data on channel refChan from the 
%               data on all OTHER channels (i.e. the data on refChan remain unreferenced themselves). default = empty
% MS 2011-02-25 compute preStimdurs, array of nstim*nrepeats, which contains the durations 
%               in between the stimuli and the stimuli prior to them

% FOUND BUG (AS): Works only with Michigan electrodes (please elaborate??)
% FOUND BUG (ND): Looks to me as if calculation of preStimdurs is wrong.
% Whoever uses this might want to check.

%% Initialization
global DIRS
global pepNEV

if ~exist('DIRS','var') || isempty(DIRS) || ~isfield(DIRS,'Cerebus')
    SetDefaultDirs;
end

if nargin < 3
    ExptTag      = varargin{1};
    animal       = ExptTag.animal;
    iseries      = ExptTag.iseries;
    iexp         = ExptTag.iexp;
else
    animal       = varargin{1};
    iseries      = varargin{2};
    iexp         = varargin{3};
end

if nargin > 3
    refChan = varargin{4};
else 
    refChan = [];
end

if nargin == 5
    preStimTime   = varargin{5}(1);
    postStimTime  = varargin{5}(2);
else
    preStimTime   = 0;
    postStimTime  = 0;
end
   
% added by ND 20101020 because the time stamps are recorded at 30000
% samples per second. used in calculation of stimdurs
nevSamplingRateInKHZ = 30;

%% See if an expt file has already been saved

AnimalDir = fullfile(DIRS.data,animal);
SeriesDir = fullfile(AnimalDir,num2str(iseries));
ExpDir    = fullfile(SeriesDir,num2str(iexp));
FileName  = fullfile(ExpDir, 'Expt_CerebusTraces.mat');
% 2011-05-31 ND made this change because exist returns 1 if it is a
% variable in the workspace and it returns a 2 if it is a file in the path
% if exist(FileName,'file')==1 && preStimTime==0% 20110511 ts added '==1'
if exist(FileName,'file')==2 && preStimTime==0% 20110511 ts added '==1'
    fprintf('Reading existing expt file... ');
    load( FileName, 'expt' );
    fprintf('done\n');
    return
end

%% Prepare to read Cerebus file

% ND 2010-10
% put in if statement so that won't break if doesn't find file. and so that
% it sequentially checks where the data should be and if they are stored by
% experiment or by repeat.

% first check to see if the data are saved by experiment and are stored by animal
byexptfname_oldway = [BuildFileName(fullfile(DIRS.Cerebus,animal), sprintf('u%03d',iseries), sprintf('%03d',iexp)) '.ns*'];
byexptfilesFound_oldway = dir(byexptfname_oldway);
% next check to see if the data are saved by experiment but in a subdirectory
byexptfname_newway = [BuildFileName(fullfile(DIRS.Cerebus,animal,num2str(iseries),num2str(iexp)),animal,iseries,iexp) '.ns*'];
byexptfilesFound_newway = dir(byexptfname_newway);
% finally check to see if the data are saved by repeat in a subdirector
byrptfname = [BuildFileName(fullfile(DIRS.Cerebus,animal,num2str(iseries),num2str(iexp)),animal,iseries,iexp) '_*.ns*'];
byrptfilesFound = dir(byrptfname);
if ~isempty(byexptfilesFound_oldway)
    byexptflag = 1;
    nfiles = 1;
    % we might be reading an .ns3 file, an .ns4 file or an .ns5 file; this gets the right extension
    byexptfname_file = byexptfilesFound_oldway.name;
    fname{1} = BuildFileName(fileparts(byexptfname_oldway),byexptfname_file);
    % if not there, check to see if the data are saved by experiment but saved in subdirectories of animal/series/expt
elseif ~isempty(byexptfilesFound_newway)
    byexptflag = 1;
    nfiles = 1;
    % we might be reading an .ns3 file, an .ns4 file or an .ns5 file; this gets the right extension
    byexptfname_file = byexptfilesFound_newway.name;
    fname{1} = BuildFileName(fileparts(byexptfname_newway),byexptfname_file);
elseif ~isempty(byrptfilesFound)
    byexptflag = 0;
    nfiles = size(byrptfilesFound,1);
    byrptfname_file = {byrptfilesFound.name};
    % potential for files to be loaded with repeats out of order when have
    % more than 9 repeats because the files found with the dir function are
    % sorted so that starts with experiment 1 then
    % 10,11,12,...,2,20,21...,3, etc. this bit sorts them properly so that
    % when data are stored in expt.data, the repeats are correct.
    expfilesorder = NaN(nfiles,1);
    for ifile = 1:nfiles
        t = regexp(byrptfname_file{ifile}, '_(\d*)\.', 'tokens');
        expfilesorder(ifile) = str2num(cell2mat(t{1})); %#ok<ST2NM>
    end
    [~,sorted_order]=sortrows([expfilesorder,(1:nfiles)']);
    fname_temp = cell(1,nfiles);
    % we might be reading an .ns3 file, an .ns4 file or an .ns5 file; this gets the right extension
    for ifile = 1:nfiles
        fname_temp{ifile} = BuildFileName(fileparts(byrptfname),byrptfname_file{sorted_order(ifile)});
    end
else
    error('Could not find any data files');
end

%% Load protocol

[p,success] = ProtocolLoad(animal,iseries,iexp);
if ~success, expt = []; return; end

% added by MC 2011-04: 
p.seqnums = p.seqnums(:,1:p.nrepeats); % removes bad endings with zeros...
% added by ND 2011-04:
if ~byexptflag
    fname = fname_temp(1:p.nrepeats); % if there is a bad ending removes from the list
    nfiles = p.nrepeats;
end

if any(p.seqnums(:) ==0),  disp('Sorry, bad sequnums'); expt = []; return; end

% alert to the possible absence of mat files
% 2010-09 ND took it out.
% probably should take this out because searching for files created by
% multispike and we don't use multispike anymore.
% ProtocolAbsentMatFiles(p);

%% Open .nev file, load photodiode data

% initialize the stimulus duration
stimdurs = zeros( p.nstim,  p.nrepeats );

for ifile = 1:nfiles
    [nevopen_outcome,foo,foo] = nevopen([fname{ifile}(1:end-2) 'ev']); %#ok<NASGU,*ASGLU>
    if nevopen_outcome < 0
        error('ExptLoadBlackRock:NevOpenError','Could not open %s.',fname{ifile});
    end

    % photodiode times that are closer than 5ms at 30k acquisition rate are
    % eliminated
    PDtimestamps  = pepNEV.sync.timestamps;
    PDtimestamps(find(diff(PDtimestamps)<5*nevSamplingRateInKHZ) + 1) = [];
    
    % photodiode times in units of clock ticks
    stimOnTimes  = PDtimestamps(1:2:end);
    stimOffTimes = PDtimestamps(2:2:end);
    
    % error checking
     if diff([length(stimOnTimes) length(stimOffTimes) max(p.seqnums(:,ifile))]) ~= 0
        error('umatched number of start/stop sync times or number of trials for %s', fname{ifile});
    end

    %% get timestamp (Greenwich mean time)
    
    timestamp = pepNEV.TimeIndex';
    % BUG HERE timestamp(2) = timestamp(3)*1+timestamp(2); % to get the month as 1 number
    timestamp = datestr(timestamp([1 2 4 5 6 7]));

    %% check number of samples for each stimulus
    % ND 20101021 changed because sometimes want to do by repeat so have to reduce the seqnums to one that is per repeat
    if byexptflag
        nsamples = stimOffTimes(p.seqnums) - stimOnTimes(p.seqnums);
    else % by repeat
        nsamples = stimOffTimes(mod(p.seqnums(:,ifile)-1,p.nstim)+1) - stimOnTimes(mod(p.seqnums(:,ifile)-1,p.nstim)+1);
    end
    nsamples = unique(nsamples);
    
    if size(nsamples,1) > 1 || size(nsamples,2) > 1
        fprintf('\nwarning: mismatch in number of samples between stimuli\n');
        fprintf('%d ', nsamples);
        fprintf('\n\n');
    end

    %% Open .ns* file and load data

    [nsopen_outcome,SamplingRateInKHZ,nchan] = nsopen(fname{ifile});
    period = nevSamplingRateInKHZ / SamplingRateInKHZ; 
    % because sync times recorded at 30KHz but the continuous data could be
    % recorded at some other resolution

    if nsopen_outcome < 0
        error('ExptLoadBlackRock:NsOpenError','Could not open %s.',fname{ifile});
    end

    % load data

    fprintf('Loading data %s...',fname{ifile});

    % waveforms = pepNEV.ns.Data.data;

    fprintf('done\n');
    
    %% Define pre- and post stimulus durations in number of samples
    sampPreStimTime  = preStimTime  * SamplingRateInKHZ * 1000; 
    sampPostStimTime = postStimTime * SamplingRateInKHZ * 1000;

    %% Get Michigan Layout
    [arrayLayout] = MichiganGetLayout(animal, iseries, nchan);

    %% fill expt.data

    fprintf('Writing continuous data to expt.data');

    for istim = 1 : p.nstim
        fprintf('\nStimulus %d', istim);
        if byexptflag
            nrepeats = 1:p.nrepeats;
        else % by repeats
            nrepeats = ifile;
        end
        for irepeat =  nrepeats % because loading by repeat anyway and already looping by irpt
            fprintf('.');
            % find the appropriate segment
            if byexptflag
                t0 = stimOnTimes(p.seqnums(istim,irepeat)) - sampPreStimTime;
                t1 = stimOffTimes(p.seqnums(istim,irepeat)) + sampPostStimTime;
                treal0 = stimOnTimes(p.seqnums(istim,irepeat));
                if p.seqnums(istim,irepeat) > 1
                    tprev1 = stimOffTimes(p.seqnums(istim,irepeat)-1);
                else
                    tprev1 = 0;
                end
            else % by repeats
                t0 = stimOnTimes(mod(p.seqnums(istim,irepeat)-1,p.nstim)+1) - sampPreStimTime;
                t1 = stimOffTimes(mod(p.seqnums(istim,irepeat)-1,p.nstim)+1) + sampPostStimTime;
                treal0 = stimOnTimes(mod(p.seqnums(istim,irepeat)-1,p.nstim)+1);
                if (mod(p.seqnums(istim,irepeat)-1,p.nstim)) > 1
                    tprev1 = stimOffTimes(mod(p.seqnums(istim,irepeat)-1,p.nstim));
                else
                    tprev1 = 0;
                end
            end
            stimdurs(istim, irepeat) = (t1-t0) /(nevSamplingRateInKHZ*1000);
            preStimdurs(istim, irepeat) = (treal0-tprev1) /(nevSamplingRateInKHZ*1000);
            
            if stimdurs(istim,irepeat) > 0.01
                % convert to units of samples
                i0 = ceil(t0/period);
                i1 = floor(t1/period);
                
                for ichan = 1 : nchan
                    if ichan<=length(arrayLayout)
                        % remap the sites according to probe layout
                        jchan = arrayLayout(ichan);
                    else
                        % don't remap a thing
                        jchan = ichan;
                    end
                    % get correct data snippet
                    expt.data{ichan}{istim,irepeat} = downsample(single(pepNEV.ns.Data.data(jchan,i0:i1)),2);
                end
            end
        end
    end
    
    %% internally reference data (optional)
    if ~isempty(refChan)
        fprintf('Internally reference data');
        for istim = 1 : p.nstim
            if byexptflag
                nrepeats = 1:p.nrepeats;
            else % by repeats
                nrepeats = ifile;
            end
            for irep = nrepeats
                for ichan = 1 : nchan 
                    if ichan ~= refChan
                        expt.data{ichan}{istim,irep} = expt.data{ichan}{istim,irep} - expt.data{refChan}{istim,irep};
                    end
                end
            end
        end
    end
        
    fprintf('\t\t\t done.\n');
    
end

SamplingRateInKHZ = SamplingRateInKHZ/2; % we downsampled

%% ok to close now

nevclose;

%% define filter coefficients for later spike sorting
switch (SamplingRateInKHZ)
    case (1) % ns3 file, originally sampled at 2000 s/sec
        filter_bounds = [100 499];
    case (5) % ns4 file, originally sampled at 10000 s/sec
        filter_bounds = [500 2300];
    case (15) % ns5 file, originally sampled at 30000 s/sec
        filter_bounds = [500 7000];
end
[b,a] = butter(3,2*filter_bounds/SamplingRateInKHZ/1000);
%Bandpass Butterworth filter between 500 Hz and 7 kHz

%% assign arbitrary threshold to data
threshvalues = zeros(nchan,1);
threshsigns  = -ones(nchan,1);

% MC 2011-Apr-07 changed threshold from -2 to -3 
for ichan = 1 : nchan
    threshvalues(ichan) = -3*std( filter(b,a,expt.data{ichan}{1,1}) );
end


%% Write expt struct for output
expt.animal           = animal;
expt.iseries          = iseries;
expt.iexp             = iexp;
expt.nrepeats         = p.nrepeats;
expt.nstim            = p.nstim;
expt.stimdurs         = stimdurs;
expt.preStimdurs      = preStimdurs;
expt.prestimtime      = preStimTime;
expt.poststimtime     = postStimTime;
expt.nchans           = nchan;
expt.refChan          = refChan;
expt.data             = expt.data';
expt.timestamp        = timestamp;
expt.channames        = (1:nchan)+500;
expt.unitspervolt     = NaN;
expt.samplerate       = SamplingRateInKHZ * 1000;
expt.threshsigns      = threshsigns';
expt.threshvalues     = threshvalues';
expt.ChanDescriptions = {};
expt.timestamps       = {};
expt.filtercoeffs     = {b,a};
expt.ecg              = [];
expt.eeg              = [];
expt.Source           = 'CerebusTraces';

%% Save the experiment file

ExptSave(expt);

