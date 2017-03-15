function [ expt, nnewrepeats ] = ExptLoadMichigan(animal, iseries, iexp, expt)
% ExptLoadMichigan loads or updates a Michigan experiment
%
% expt = ExptLoadMichigan(ExptTag ) loads the expt structure associated with
% ExptTag.animal, ExptTag.iseries, ExptTag.iexp.
%
% expt = ExptLoadMichigan(animal, iseries, iexp) lets you specify those entries
%
% expt = ExptLoadMichigan(animal, iseries, iexp, expt) updates with the repeats (if
% any) that were not already loaded 
%
% expt = ExptLoadMichigan(animal, iseries, iexp, batchflag), lets you specify that
% you are working in batch (so, no popup windows) batchflag:=
% {'batch','quiet'}
%
% [ expt, nnewrepeats ] = ... tells you how many new repeats were loaded.
%
% 2000-01 Matteo Carandini 2000-09 MC only mat data 
% 2001-11 MC merged with ExptUpdate 
% 2003-04 MC changed behavior when attempting to load a bad datafile 
% 2003-08 VB added 'batch' and 'quiet' flags. 
% 2004-03 VM corrected eeg bug 
% 2007-06 MC added support for one single argument 
% 2007-09 MC added hack to load Michigan files -- ignores multispike data!!!
% 2007-09 MC split into ExptLoadMichigan
% 2007-11 MC added the "try-catch" when loading to avoid mem overflow
% 2009-05 MC added call to ProtocolAbsentMatFiles
% 2009-11 MC tried to make it resist the absence of mat files

if nargin == 1
    ExptTag = animal;
    animal = ExptTag.animal;
    iseries = ExptTag.iseries;
    iexp = ExptTag.iexp;
end

global DIRS

[protocol, success] = ProtocolLoad(animal, iseries, iexp);
if ~success, expt = []; return; end

% alert to the possible absence of mat files
ProtocolAbsentMatFiles(protocol);

verbose = 2;
if nargin >=4 && ischar(expt)
    if strcmp(expt,'quiet') % as quiet as possible
        verbose =0;
    elseif strcmp(expt,'batch') 
        % dialog boxes directed to stdout
        verbose = 1;
    end
end

if nargin < 4 || ischar(expt)
   expt = ExptInitialize;  
   expt.animal      = animal;  
   expt.iseries     = iseries;  
   expt.iexp        = iexp;
   expt.nstim       = protocol.nstim;
   expt.nrepeats    = 0;
end

oldnrepeats = expt.nrepeats;

% ----------------------- load the data
if verbose>0;fprintf('Reading MAT data (%d stimuli)', protocol.nstim);end;

flgError = 0;
flgTrodes = 0; % assume there were no trode traces
flgCamera = 0; % assume there was no camera

% just an initial allocation:
dd          = cell(protocol.nrepeats, protocol.nstim);
frametoggle	= cell(protocol.nrepeats, protocol.nstim);

irepeat 	= oldnrepeats;
nrepeats	= oldnrepeats;

while nrepeats == irepeat
    irepeat = irepeat+1; % let's try to read the next repeat
    
    MemInfo = whos('dd');
    if MemInfo.bytes > 2^30
        disp('Exceeded 1 GB of data -- not loading any further repeats');
        flgError = 1; %#ok<NASGU>
        break;
    end
    
    for istim = 1:protocol.nstim
               
        filename = sprintf('%s_%d_%d_%d_%d.mat',animal, iseries, iexp, irepeat, istim);
        datafile = fullfile(DIRS.data,animal,int2str(iseries),int2str(iexp),filename);
        if exist(datafile,'file')~=2
            MultispikeData = struct();
        else
            load(datafile); % we are loading stimresps
            MultispikeData = stimresps;
            clear stimresps
        end
        
        if isfield(MultispikeData,'frametoggle')
            flgCamera = 1;
            frametoggle{irepeat,istim} = MultispikeData.frametoggle;
        end
        
        if isfield(MultispikeData,'timestamp')
            expt.timestamps{istim,irepeat} = MultispikeData.timestamp;
        end
        
        if isfield(MultispikeData,'ecg')
            expt.ecg.traces{istim,irepeat} = double(MultispikeData.ecg.data);
            expt.ecg.rate = MultispikeData.ecg.rate;
        end
        if isfield(MultispikeData,'eeg')
            expt.eeg.traces{istim,irepeat} = double(MultispikeData.eeg.data);
            expt.eeg.rate = MultispikeData.eeg.rate;
        end

        % read the Michigan file
        MichiganFileName = sprintf('%s_%d_%d_%d_%d-Michigan.mat',...
            animal, iseries, iexp, irepeat, istim);
        MichiganDataFile = fullfile(DIRS.data,animal,int2str(iseries),int2str(iexp),MichiganFileName);
        if ~exist(MichiganDataFile,'file')
            flgError = 1;
        else
            flgTrodes = 1; 
            try
                load(MichiganDataFile,'stimresps'); 
            catch %#ok<CTCH>
                % probably went out of memory
                disp('Probably exceeded available memory -- not loading any further repeats');
                flgError = 1;
                break;
            end
            MichiganData = stimresps; clear stimresps; % rename
            
            dd{irepeat,istim} = MichiganData.data;
            MichiganData.threshsigns = -1*ones(16,1);
            MichiganData.threshvalues = 0.85*double(min(dd{irepeat,istim},[],2));
            MichiganData.channames = 101:116;
            MichiganData.MaxV = 1;
            % MichiganData.filtercoeffs = {a,b}; % is this the right format??
%             NumBytes = numel(MichiganData.data)*2;
%             if NumBytes > feature('memstats') % loading one more will lead to out of memory 
%                 disp('Exceeded available memory -- not loading any further repeats');
%                 flgError = 1;
%                 break;
%             end
        end
                
    end
    if flgError
        break;
    end
    nrepeats = irepeat; % we succeeded in reading this repeat
    if verbose>0, fprintf('.');end;
end

% the following code is complicated because expt.timestamps can be empty
if isfield(expt,'timestamps')
    dates = expt.timestamps(1:end,1); % the times of the first stimulus
    for idate = 1:length(dates)
        if ~isempty(dates{idate})
            dates{idate} = datenum(dates{idate});
        end
    end
    expt.timestamp = datestr(min([dates{:}]));
end

nnewrepeats = nrepeats-oldnrepeats;
if verbose>0, fprintf(' done (%d / %d repeats)\n',nnewrepeats,nrepeats);end;

%%

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

%%
% this may not work in those files in which we only have frametoggle and no
% electrodes.

expt.samplerate     = MichiganData.samplerate;
expt.unitspervolt	= NaN;

if flgTrodes
    
    % TRODE DATA WERE ACQUIRED
    
    expt.threshsigns 	= MichiganData.threshsigns;
    expt.threshvalues	= MichiganData.threshvalues;
    expt.channames 	    = MichiganData.channames;
    if isfield(MichiganData,'filtercoeffs')
        expt.filtercoeffs = MichiganData.filtercoeffs;
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
                expt.data{ichan}{istim,irepeat} = dd{irepeat,istim}(ichan,:); 
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
    if isnan(expt.channames)
        disp('Actual channel names were not saved. Setting them arbitrarily to 1, 2,...');
        expt.channames = 1:expt.nchans;
    end

end

if flgCamera
    
    % CAMERA DATA WERE ACQUIRED
        
    for irepeat = (expt.nrepeats+1):nrepeats
        for istim = 1:protocol.nstim
            expt.frametoggle{istim,irepeat} = frametoggle{irepeat,istim};
            % expt.stimdurs(istim,irepeat) = length(frametoggle{irepeat,istim})/MultispikeData.samplerate;
        end
    end
    
end

expt.nrepeats = nrepeats; 



