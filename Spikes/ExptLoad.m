function [ expt, nnewrepeats ] = ExptLoad(animal, iseries, iexp, expt)
% ExptLoad loads or updates an experiment
%
% expt = ExptLoad(ExptTag ) loads the expt structure associated with
% ExptTag.animal, ExptTag.iseries, ExptTag.iexp.
%
% expt = ExptLoad(animal, iseries, iexp) lets you specify those entries
%
% expt = ExptLoad(animal, iseries, iexp, expt) updates with the repeats (if
% any) that were not already loaded 
%
% expt = ExptLoad(animal, iseries, iexp, batchflag), lets you specify that
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
% 2007-09 MC reverted to older stuff and created ExptLoadMichigan
% 2007-10 MC frametoggle is now uint16
% 2009-05 MC added call to ProtocolAbsentMatFiles
% 2010-06 MC deals with absent mat files

% global DIRS

%% Parse the parameters

if nargin == 1
    ExptTag = animal;
    animal = ExptTag.animal;
    iseries = ExptTag.iseries;
    iexp = ExptTag.iexp;
end

verbose = 2;
if nargin >=4 && ischar(expt)
    switch expt
        case 'quiet'
            verbose =0;
        case 'batch'
            verbose = 1; % dialog boxes directed to stdout
    end
    expt = [];
end

if nargin < 4
   expt = [];
end

%% load basic info

[protocol, success] = ProtocolLoad(animal, iseries, iexp);
if ~success, expt = []; return; end

if isempty(expt)
   expt = ExptInitialize;  
   expt.animal      = animal;  
   expt.iseries     = iseries;  
   expt.iexp        = iexp;
   expt.nstim       = protocol.nstim;
   expt.nrepeats    = 0;
end

oldnrepeats = expt.nrepeats;

% alert to the possible absence of mat files
[MatFileNames, AbsentFiles] = ProtocolMatFiles(protocol);

nrepeats = size(MatFileNames,2);

%% Load the data

if verbose>0; fprintf('Reading MAT data (%d stimuli)', protocol.nstim); end;

flgTrodes = 0; % assume there were no trode traces
flgCamera = 0; % assume there was no camera

dd          = cell(protocol.nstim,nrepeats);
frametoggle	= cell(protocol.nstim,nrepeats);

for irepeat = (oldnrepeats+1):nrepeats
    
    % let's try to read this repeat
    
    if all(AbsentFiles(:,irepeat))
        fprintf(1,'\n WARNING: No data files for entire repeat %d. Stopping the loading.\n', irepeat);
        nrepeats = irepeat;
        break;
    end
    
    for istim = 1:protocol.nstim
        
        datafile = MatFileNames{istim, irepeat};

        if isempty(datafile)
            fprintf(1,'\n WARNING: No data file for stim %d repeat %d\n', istim, irepeat);
            
        else
            load(datafile); % we are loading stimresps
            
            if isfield(stimresps,'data') && ~isempty(stimresps.data)
                flgTrodes = 1;
                dd{istim,irepeat} = stimresps.data;
            end
             
            if isfield(stimresps,'frametoggle')
                flgCamera = 1;
                frametoggle{istim,irepeat} = uint16(stimresps.frametoggle);
            end
            
            if isfield(stimresps,'ecg')
                expt.ecg.traces{istim,irepeat} = double(stimresps.ecg.data);
                expt.ecg.rate = stimresps.ecg.rate;
            end
            
            if isfield(stimresps,'eeg')
                expt.eeg.traces{istim,irepeat} = double(stimresps.eeg.data);
                expt.eeg.rate = stimresps.eeg.rate;
            end
                
            expt.timestamps{istim,irepeat} = stimresps.timestamp;
        end
        
    end
   
    if verbose>0, fprintf('.'); end;
end

% the following code is complicated because expt.timestamps can be empty
if isfield(expt,'timestamps');
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

%%

% this may not work in those files in which we only have frametoggle and no
% electrodes.

expt.samplerate     = stimresps.samplerate;
expt.unitspervolt	= 32768/stimresps.MaxV;

expt.AbsentFiles = AbsentFiles;

if flgTrodes
    
    % TRODE DATA WERE ACQUIRED
    
    expt.threshsigns 	= stimresps.threshsigns;
    expt.threshvalues	= stimresps.threshvalues;
    expt.channames 	    = stimresps.channames;
    if isfield(stimresps,'filtercoeffs')
        expt.filtercoeffs = stimresps.filtercoeffs;
    else
        expt.filtercoeffs = [];
    end
    
    %---------------------------------------------------------
    
    if expt.nrepeats == 0 % we are starting from scratch
        expt.nchans = size(dd{1,nrepeats},1); 
        expt.data = cell(expt.nchans,1); 
    end
     
    if expt.nchans == 0
        expt = [];
        return
    end
    
    for irepeat = (expt.nrepeats+1):nrepeats
        
        for istim = 1:protocol.nstim
            
            if isempty(dd{istim,irepeat})
                for ichan = 1:expt.nchans
                    expt.data{ichan}{istim,irepeat} = [];
                end
                expt.stimdurs(istim,irepeat) = 0;
            
            else
                
                % turn the data around
                for ichan = 1:expt.nchans
                    expt.data{ichan}{istim,irepeat} = single(dd{istim,irepeat}(ichan,:))/expt.unitspervolt;
                end
                dd{istim,irepeat} = []; % to clear memory
                expt.stimdurs(istim,irepeat) = length(expt.data{ichan}{istim,irepeat})/expt.samplerate;
            end
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



