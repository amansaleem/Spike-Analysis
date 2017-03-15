function [expt,lfp] = ExptLoadCerebusTraces(varargin)
% ExptLoadCerebusTraces loads Cerebus data from an ns* file
% This function is meant for continuous data (i.e. ns* files, not nev
% files!). It unmaps the channels for Michigan probes, so for example
% channel 3 is the third from the top in a probe, and so on. 
%
% Data are downsampled by 2. So if you acquired at 30000 Hz, data are at
% 15000 Hz. Also, LFPs are always computed at 400 Hz and stored in a
% separate file.
%
% [expt,lfp] = ExptLoadCerebusTraces(ExptTag) loads data for the specified experiment
% (ExptTag should have fields animal, iseries, iexp).
%
% [expt,lfp] = ExptLoadCerebusTraces(animal,iseries,iexp) does the same.
%
% [expt,lfp] = ExptLoadCerebusTraces(animal,iseries,iexp,refChan) lets you specify a
% channel to be used as internal reference. Channel refChan itself is left
% untouched. This channel is the one after the probe has been unmapped,
% i.e., 3 means to reference from channel 3 - the third from the top in an
% A1x16 michigan probe. Default: [].
%
% [expt,lfp] = ExptLoadCerebusTraces(animal,iseries,iexp,..., [pre post]) lets you
% specify the number of seconds to prepend and append to the traces. For
% example, [3 1] means 3 s of activity before stimulus onset, and 1 s of
% activity after stimulus offset. Default: [], which means [0 0].
%
% [expt,lfp] = ExptLoadCerebusTraces(ExptTag,[], [], [], ,LoadWhat) lets
% you specify what files to load/save: LoadWhat = 'Both' (DEFAULT) or 'LFPonly'.
%
% NOTICE THAT TRACES ARE DOWNSAMPLED BY A FACTOR OF 2 UPON LOADING

% MS 2010-03-22 created (adapted from UtahGetLFP)
% MS 2010-03-24 minor changes to calculating sample rate and stimulus durations
% MC 2010-06-14 changed the order of opening ns5 and nev (trouble with memmap)
% MC 2010-06-14 introduced a downsample by a factor of 2
% MC 2010-06-16 attempted to fix bug in timestamp (search for BUG)
% MC 2010-06-16 saves the data on zserver as a matlab file (an "expt")
% MC 2010-06-16 From then on, when you call it, it loads the mat file directly. 
% AS 2010-07-19 Added the functionality to load a Michigan map file
% ND 2010-10-10 took out call to ProtocolAbsentMatfiles (we no longer use Multispike)
% ND 2010-10-21 Cerebus files named more appropriately and saved by repeat rather than by experiment.
% MS 2010-10-25 added another input argument: StimTime [pre post]
% ND 2010-12-09 Fixed bug by AS 2010-07-19 (last site was not remapped)
% ND 2010-12-09 Fixed bug by ND 2010-10-21 (weird assignment of repeats).
% MS 2010-12-14 when StimTime exists (see 2010-10-25), this code is always run, even when an expt file already exists 
% AA 2011-02-01 deleted photodiode time stamps (PDtimestamps) that are closer than 5 ms
% MS 2011-02-23 added another input argument: refChan, that subtracts the data on channel refChan 
% MS 2011-02-25 compute preStimdurs, array of nstim*nrepeats, which contains the pre-stim durations 
% ND 2011-07-07 creates an exptfp structure that is a lot like the expt but for lfps
% ND 2011-07-07 created optional input argument of samplerate for FPs
% MC 2011-07-12 simplified and streamlined 
% MS 2011-07-18 fixed small bug in computing preStimdurs when data saved by repeat
% AS 2011-07-21 de-bugged loading existing files
% ND 2012-03-28 added optional input argument to only load FP data

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
% refChan = input('ref Channel number or skip (Enter)>>');%arrayLayout is considered

if nargin > 4
    PrePost = varargin{5};
else
    PrePost = [];
end

if isempty(PrePost)
    preStimTime   = 0;
    postStimTime  = 0;
else
    preStimTime   = PrePost(1);
    postStimTime  = PrePost(2);
end

LoadWhat = 'Both'; % default
if nargin > 5
    LoadWhat = varargin{6};
    if ~any(strcmp({'Both','LFPonly'},LoadWhat))
        error('Argument LoadWhat can only be ''Both'' or ''LFPonly''');
    end
end

% time stamps are recorded at 30000 Hz in the NEV file
nevSamplingRateInKHZ = 30;

desired_lfp_samplerate = 400; % Hz. 
 
%% See if an expt file has already been saved

AnimalDir  = fullfile(DIRS.data,animal);
SeriesDir  = fullfile(AnimalDir,num2str(iseries));
ExpDir     = fullfile(SeriesDir,num2str(iexp));

FileName    = fullfile(ExpDir, 'Expt_CerebusTraces.mat');
FileNameLFP = fullfile(ExpDir, 'Expt_CerebusTraces_LFP.mat' );

if exist(FileName,'file')==2 && preStimTime==0 && exist(FileNameLFP,'file')==2
    
    expt = [];
    if strcmp(LoadWhat,'Both')
        fprintf('Reading existing expt file... ');
        load( FileName, 'expt' );
        fprintf('done\n');
    end
    
    fprintf('Reading existing lfp file... ');
    load( FileNameLFP, 'lfp' );
    fprintf('done\n');
return % AS 2011-07-21, re-loads the data from the ns5 file without this statement
end

%% Prepare to read Cerebus file

DataType = ''; % can be 'ByExptInAnimalDir', 'ByExptInExptDir', 'ByRepeat'

% first check to see if the data are saved by experiment and are stored by animal
FilePattern = [BuildFileName(fullfile(DIRS.Cerebus,animal), sprintf('u%03d',iseries), sprintf('%03d',iexp)) '.ns*'];
if ~isempty(dir(FilePattern)), DataType = 'ByExptInAnimalDir'; end

% next check to see if the data are saved by experiment but in a subdirectory
FilePattern = [BuildFileName(fullfile(DIRS.Cerebus,animal,num2str(iseries),num2str(iexp)),animal,iseries,iexp) '.ns*'];
if ~isempty(dir(FilePattern)), DataType = 'ByExptInExptDir'; end

% finally check to see if the data are saved by repeat
FilePattern = [BuildFileName(fullfile(DIRS.Cerebus,animal,num2str(iseries),num2str(iexp)),animal,iseries,iexp) '_*.ns*'];
if ~isempty(dir(FilePattern)), DataType = 'ByRepeat'; end

%%

switch DataType
    
    case 'ByExptInAnimalDir'
        FilePattern = [BuildFileName(fullfile(DIRS.Cerebus,animal), sprintf('u%03d',iseries), sprintf('%03d',iexp)) '.ns*'];
        FoundFiles = dir(FilePattern);
        nfiles = 1;
        FullFileName{1} = BuildFileName(fileparts(FilePattern),FoundFiles.name);
        
    case 'ByExptInExptDir'
        FilePattern = [BuildFileName(fullfile(DIRS.Cerebus,animal,num2str(iseries),num2str(iexp)),animal,iseries,iexp) '.ns*'];
        FoundFiles = dir(FilePattern);
        nfiles = 1;
        FullFileName{1} = BuildFileName(fileparts(FilePattern),FoundFiles.name);
        
    case 'ByRepeat'
        FilePattern = [BuildFileName(fullfile(DIRS.Cerebus,animal,num2str(iseries),num2str(iexp)),animal,iseries,iexp) '_*.ns*'];
        FoundFiles = dir(FilePattern);
        nfiles = size(FoundFiles,1);
        byrptfname_file = {FoundFiles.name};
        
        % sort files so repeat 9 is before repeat 10
        expfilesorder = NaN(nfiles,1);
        for ifile = 1:nfiles
            t = regexp(byrptfname_file{ifile}, '_(\d*)\.', 'tokens');
            expfilesorder(ifile) = str2num(cell2mat(t{1})); %#ok<ST2NM>
        end
        [~,sorted_order]=sortrows([expfilesorder,(1:nfiles)']);
        
        FullFileName = cell(1,nfiles);
        for ifile = 1:nfiles
            FullFileName{ifile} = BuildFileName(fileparts(FilePattern),byrptfname_file{sorted_order(ifile)});
        end
        
    otherwise
        error('Could not find any data files');
end

%% Load protocol

[p,success] = ProtocolLoad(animal,iseries,iexp);
if ~success, expt = []; return; end

% added by MC 2011-04:
p.seqnums = p.seqnums(:,1:p.nrepeats); % removes bad endings with zeros...
% added by ND 2011-04:
if strcmp(DataType,'ByRepeat')
    FullFileName = FullFileName(1:p.nrepeats); % if there is a bad ending removes from the list
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
lfp.data = cell(p.nstim,p.nrepeats);

 
for ifile = 1:nfiles
    [nevopen_outcome,foo,foo] = nevopen([FullFileName{ifile}(1:end-2) 'ev']); %#ok<NASGU,*ASGLU>
    if nevopen_outcome < 0
        error('ExptLoadBlackRock:NevOpenError','Could not open %s.',FullFileName{ifile});
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
        error('umatched number of start/stop sync times or number of trials for %s', FullFileName{ifile});
    end
    
    %% get timestamp (Greenwich mean time)
    
    timestamp = pepNEV.TimeIndex';
    % BUG HERE timestamp(2) = timestamp(3)*1+timestamp(2); % to get the month as 1 number
    timestamp = datestr(timestamp([1 2 4 5 6 7]));
    
    %% check number of samples for each stimulus

    % MC 2011-07 commented this out. In general different number of samples
    % should not be a problem. 
    
%     if strcmp(DataType,'ByRepeat')
%         nsamples = stimOffTimes(mod(p.seqnums(:,ifile)-1,p.nstim)+1) - stimOnTimes(mod(p.seqnums(:,ifile)-1,p.nstim)+1);
%     else
%         nsamples = stimOffTimes(p.seqnums) - stimOnTimes(p.seqnums);
%     end
%     nsamples = unique(nsamples);
%     
%     if size(nsamples,1) > 1 || size(nsamples,2) > 1
%         fprintf('\nwarning: mismatch in number of samples between stimuli\n');
%         fprintf('%d ', nsamples);
%         fprintf('\n\n');
%     end
    
    %% Open .ns* file and load data
    
    [nsopen_outcome,SamplingRateInKHZ,nchan] = nsopen(FullFileName{ifile});
    period = nevSamplingRateInKHZ / SamplingRateInKHZ;
    % because sync times recorded at 30KHz but the continuous data could be
    % recorded at some other resolution
    
    if nsopen_outcome < 0
        error('ExptLoadBlackRock:NsOpenError','Could not open %s.',FullFileName{ifile});
    end
    
    % load data
    
%     fprintf('Loading data %s...',FullFileName{ifile});
%     
%     % waveforms = pepNEV.ns.Data.data;
%     
%     fprintf('done\n');
    
    %% Define pre- and post stimulus durations in number of samples
    sampPreStimTime  = preStimTime  * SamplingRateInKHZ * 1000;
    sampPostStimTime = postStimTime * SamplingRateInKHZ * 1000;
    
    %% Get Michigan Layout
    [arrayLayout] = MichiganGetLayout(animal, iseries, nchan);
    
    %% get ready for LFPs
    
    % dowsample from original samplerate to desired_lfp_samplerate
    downsample_factor = floor(15000/desired_lfp_samplerate);
    lfp.samplerate = 15000/downsample_factor;

   %% fill expt.data
    
%    fprintf('Writing continuous data to expt.data');
    
    for istim = 1 : p.nstim
        fprintf('\nStimulus %d', istim);
        if strcmp(DataType,'ByRepeat')
            RepeatList = ifile;
        else
            RepeatList = 1:p.nrepeats;
        end
        for irepeat =  RepeatList 
            fprintf(' Repeat %d ', irepeat );
            % find the appropriate segment
            if strcmp(DataType,'ByRepeat')
                t0 = stimOnTimes(mod(p.seqnums(istim,irepeat)-1,p.nstim)+1) - sampPreStimTime;
                t1 = stimOffTimes(mod(p.seqnums(istim,irepeat)-1,p.nstim)+1) + sampPostStimTime;
                treal0 = stimOnTimes(mod(p.seqnums(istim,irepeat)-1,p.nstim)+1);
                if (mod(p.seqnums(istim,irepeat)-1,p.nstim)) > 0
                    tprev1 = stimOffTimes(mod(p.seqnums(istim,irepeat)-1,p.nstim));
                else
                    tprev1 = 0;
                end
           
            else
                t0 = stimOnTimes(p.seqnums(istim,irepeat)) - sampPreStimTime;
                t1 = stimOffTimes(p.seqnums(istim,irepeat)) + sampPostStimTime;
                treal0 = stimOnTimes(p.seqnums(istim,irepeat));
                if p.seqnums(istim,irepeat) > 1
                    tprev1 = stimOffTimes(p.seqnums(istim,irepeat)-1);
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
                    
                    lfp.data{istim,irepeat}(ichan,:) = resample( ...
                        double(expt.data{ichan}{istim,irepeat}),1,downsample_factor);
                    
                end
            end
        end
    end
    
    %% internally reference data (optional)
    if ~isempty(refChan)
        fprintf('Internally reference data');
        for istim = 1 : p.nstim
            if strcmp(DataType,'ByRepeat')
                RepeatList = ifile;
            else
                RepeatList = 1:p.nrepeats;
            end
            for irep = RepeatList
                for ichan = 1 : nchan
                    if ichan ~= refChan
                        expt.data{ichan}{istim,irep} = expt.data{ichan}{istim,irep} - expt.data{refChan}{istim,irep};
                    end
                end
            end
        end
    end
    
    fprintf('\n\n');
    
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

%% Write expt struct for output
lfp.animal           = animal;
lfp.iseries          = iseries;
lfp.iexp             = iexp;
lfp.nrepeats         = p.nrepeats;
lfp.nstim            = p.nstim;
lfp.stimdurs         = stimdurs;
lfp.preStimdurs      = preStimdurs;
lfp.prestimtime      = preStimTime;
lfp.poststimtime     = postStimTime;
lfp.nchans           = nchan;
lfp.channames        = (1:nchan)+500;
lfp.ChanDescriptions = {};
lfp.timestamps       = {};
lfp.Source           = 'CerebusTraces';

%% Save the experiment file together with the LFP data

ExptSave(expt,lfp);
