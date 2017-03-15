function expt = ExptLoadUtah(varargin)
% EXPTLOADUTAH.M: takes an ExptTag and a list of channels (eg [6 26])
% and returns an "expt" structure.
% 
% Usage:
% The last argument is a vector with containing channel numbers
% Can pass an ExptTag (has fields animal, series, iexp) as the 1st argument:
% expt = ExptLoadUtah(ExptTag,chanList)
% 
% Can also as the 1st 3 arguments the animal, series, iexp values themselves:
% expt = ExptLoadUtah(animal,series,iexp,chanList)
% 
% AZ 2010-02-05 Created

%% Initialization
global DIRS

verbose = true;

if ~exist('DIRS','var') || isempty(DIRS) || ~isfield(DIRS,'Cerebus')
   SetDefaultDirs;
end

if     nargin == 2
   ExptTag      = varargin{1};
   chanList     = varargin{2};
   animal       = ExptTag.animal;
   iseries      = ExptTag.iseries;
   iexp         = ExptTag.iexp;
elseif nargin  > 2
   animal       = varargin{1};
   iseries      = varargin{2};
   iexp         = varargin{3};
    if nargin  > 3
       chanList = varargin{4};
    else
       error('Please provide a chanList.')
    end
end
% Transpose for loops later on, if necessary
if size(chanList,1) > size(chanList,2)
   chanList = chanList';
end

% set fname to filename of .nev file we want to open
fname = sprintf('%s%su%03i_%03i.nev',[DIRS.Cerebus filesep],...
   [animal filesep],iseries,iexp);

% choose threshold (32 is a good default)
threshold = 0;

timestamps    = cell(numel(chanList),1);
Waveforms_raw = timestamps;
Wav_realigned = timestamps;
iXraw         = timestamps;

%% Load protocol
[protocol, success] = ProtocolLoad(animal,iseries,iexp);
if ~success, expt = []; return; end

% alert to the possible absence of mat files
ProtocolAbsentMatFiles(protocol);

%% Open .nev file, load data
[nevopen_outcome,SamplingRateInKHZ,nChans] = nevopen(fname);
if nevopen_outcome < 0
   error('ExptLoadUtah:NevOpenError','Could not open %s.',fname);
end

% LOAD DATA
fprintf('Loading data');
for elec = chanList
   fprintf('.');
   % choose elec
   [timestamps{elec}, Waveforms_raw{elec}] = nevwaves(elec);
end
fprintf(' done.\n');
nevclose;

% Rethreshold, realign
for elec = chanList
   % Wav_realigned is Waveforms_raw(:,~iXraw) combined with realigned 
   % Waveforms_raw(:,iXraw)
   % Wav_realigned(:,iXraw) are all the traces (realigned) above threshold
   % Waveforms_raw(:,iXraw) are all the raw traces (not realigned) above threshold
   [Wav_realigned{elec},iXraw{elec}] = ...
      changeThreshold(threshold,Waveforms_raw{elec},[],nevopen_outcome);
end

%% Figure out how to arrange spike data

%% MC CODE
oldnrepeats = 0;
% ----------------------- load the data
if verbose>0;fprintf('Reading MAT data (%d stimuli)', protocol.nstim);end;

flgError  = 0;
flgTrodes = 0; % assume there were no trode traces
flgCamera = 0; % assume there was no camera

dd = {};
irepeat 	= oldnrepeats;
nrepeats	= oldnrepeats;
while nrepeats == irepeat
    irepeat = irepeat+1; % let's try to read the next repeat
    for istim = 1:protocol.nstim
        filename = sprintf('%s_%d_%d_%d_%d.mat',animal, iseries, iexp, irepeat, istim);
        datafile = fullfile(DIRS.data,animal,int2str(iseries),int2str(iexp),filename);
        if exist(datafile,'file')~=2
            fprintf(' (no %s)', filename);
            flgError = 1;
            break
        end
        
        try
            load(datafile); % we are loading stimresps
        catch
            if verbose
                fprintf(1,'\n-----------------> Problems loading file %s: Aborting (there might be unread repeats)\n', datafile);
            end
            flgError = 1;
            break
        end
        
        if isfield(stimresps,'data')
            flgTrodes = 1;
            dd{irepeat,istim} = stimresps.data;
            
            % 2010-03 MC added these lines
            if isempty(dd{irepeat,istim})
                flgError = 1;
                fprintf(1,'\n WARNING: Data for stim %d repeat %d is empty!!! Further repeats will be ignored\n', istim, irepeat);
            end
            
        end
        if isfield(stimresps,'frametoggle')
            flgCamera = 1;
            frametoggle{istim,irepeat} = uint16(stimresps.frametoggle);
        end
        
        expt.timestamps{istim,irepeat} = stimresps.timestamp;
        
        if isfield(stimresps,'ecg')
            expt.ecg.traces{istim,irepeat} = double(stimresps.ecg.data);
            expt.ecg.rate = stimresps.ecg.rate;
        end
        if isfield(stimresps,'eeg')
            expt.eeg.traces{istim,irepeat} = double(stimresps.eeg.data);
            expt.eeg.rate = stimresps.eeg.rate;
        end
                
    end
    if flgError
        break;
    end
    nrepeats = irepeat; % we succeeded in reading this repeat
    if verbose>0, fprintf('.');end;
end

% the following code is complicated because expt.timestamps can be empty
if isfield(expt,'timestamp');
    dates = expt.timestamps(1:end,1);
    for idate = 1:length(dates)
        if ~isempty(dates{idate})
            dates{idate} = datenum(dates{idate});
        end
    end
    expt.timestamp = datestr(min([dates{:}]));
end

nnewrepeats = nrepeats-oldnrepeats;
if verbose>0, fprintf(' done (%d / %d repeats)\n',nnewrepeats,nrepeats);end;

%---------------------------------------------------------

if nrepeats == oldnrepeats
    % we did not load any new repeat
    if nrepeats==0
        msg = sprintf('Could not load a full repeat for Animal %s, series %d, experiment %d (interrupted? in progress?)\n',animal,iseries,iexp);
        expt = [];
    else
        msg = sprintf('Could not update with a full repeat Animal %s, series %d, experiment %d (interrupted? in progress?)\n',animal,iseries,iexp);
    end
    switch verbose
        case 0
            % do nothing
        case 1
            disp(msg);
        otherwise
            errordlg(msg,'SpikeSorter','modal');
    end
    return % leave this function
end

%---------------------------------------------------------

% this may not work in those files in which we only have frametoggle and no
% electrodes.

expt.samplerate     = stimresps.samplerate;
expt.unitspervolt	= 32768/stimresps.MaxV;

if flgTrodes
    
    % TRODE DATA WERE ACQUIRED
    
    expt.threshsigns 	= stimresps.threshsigns;
    expt.threshvalues	= stimresps.threshvalues;
    expt.channames 	   = stimresps.channames;
    if isfield(stimresps,'filtercoeffs')
        expt.filtercoeffs = stimresps.filtercoeffs;
    else
        expt.filtercoeffs = [];
    end
    
    %---------------------------------------------------------
    
    if expt.nrepeats == 0 % we are starting from scratch
       expt.nchans = size(dd{nrepeats,1},1); 
       expt.data = cell(expt.nchans,1); 
    end
    
    for irepeat = (expt.nrepeats+1):nrepeats
        
        for istim = 1:protocol.nstim
            % turn the data around  
            for ichan = 1:expt.nchans
                expt.data{ichan}{istim,irepeat} = single(dd{irepeat,istim}(ichan,:))/expt.unitspervolt;
            end
            dd{irepeat,istim} = []; % to clear memory
            expt.stimdurs(istim,irepeat) = length(expt.data{ichan}{istim,irepeat})/expt.samplerate;
        end
    end
    
    % FOR THOSE EXPTS (e.g. for early guinea pig expts) where channames was a cell array
    if iscell(expt.channames)
        expt.channames = [expt.channames{:}];
    end
    
    % HACK FOR THOSE EXPTS (hopefully only up to cat001) WHERE THE CHAN NAMES WERE NOT SAVED
    if isnumeric(expt.channames) && all(isnan(expt.channames))
        disp('Actual channel names were not saved. Setting them arbitrarily to 1, 2,...');
        expt.channames = 1:expt.nchans;
    end

end

if flgCamera   % CAMERA DATA WERE ACQUIRED
    for irepeat = (expt.nrepeats+1):nrepeats
        for istim = 1:protocol.nstim
            % expt.frametoggle{istim,irepeat} = single(frametoggle{istim,irepeat})/expt.unitspervolt;
            expt.frametoggle{istim,irepeat} = frametoggle{istim,irepeat}; % stays in uint16
            expt.stimdurs(istim,irepeat) = length(frametoggle{istim,irepeat})/expt.samplerate;
            frametoggle{istim,irepeat} = []; % important to clear memory
        end
    end
end

expt.nrepeats = nrepeats; 

%% Write expt struct for output
expt                  = ExptInitialize;

expt.animal           = animal;
expt.iseries          = iseries;
expt.iexp             = iexp;

expt.nrepeats         = 0;
expt.nstim            = NaN;
expt.stimdurs         = [];
expt.nchans           = nChans;
expt.data             = [];
expt.timestamp        = NaN;
expt.channames        = [];
expt.unitspervolt     = NaN;
expt.samplerate       = SamplingRateInKHZ*1000;
expt.threshsigns      = sign(threshold);
expt.threshvalues     =      threshold ;
expt.ChanDescriptions = {};

end