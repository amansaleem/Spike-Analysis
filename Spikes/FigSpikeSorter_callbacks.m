function FigSpikeSorter_callbacks (strarg)
% Switchyard for the FigSpikeSorter figure
%
% 2000 Matteo Carandini
% 2007-09 MC made it deal with two data types: Multispike and Michigan
% 2009-01 SK display working directory before loading discrimination parameters
% 2010-03 MC bandpass filter data acquired at high sample rates (newer data)
% 2010-03 MC folded subfunction initialize_chans into subfunction LoadedExpt
% 2010-06 MC cosmetic changes + allow saving without looking for >16 chans
% 2011-02 AS added case 'LoadContinuousData' to load VR/Cont data types

global PICK DIRS;


MatteoDebug = false;

if nargin<1
    error('FigSpikeSorter must have one argument');
end

%% If this is the first call, do nothing (PICK does not even exist yet)

if strcmp(strarg,'create')
    return
end

%% figure out which channel we are dealing with
myobj = findobj(PICK.spikesorter,'Tag','popChannel');
if ~isempty(myobj), ichan = get(myobj,'value'); end

myobj = findobj(PICK.spikesorter,'Tag','popRepeat');
if ~isempty(myobj), irepeat = get(myobj,'value'); end

%% figure out animal, series, and experiment

stranimal 	= get(findobj(PICK.spikesorter,'Tag','txtAnimal'),'String');
if iscell(stranimal), stranimal = stranimal{1}; end % there can be 2 or more txtAnimal !!!
strseries 	= get(findobj(PICK.spikesorter,'Tag','txtSeries'),'String');
strexp 		= get(findobj(PICK.spikesorter,'Tag','txtExperiment'),'String');

if ~ischar(strseries) || ~ischar(strexp)
    error('cannot read series or experiment from the spike sorter figure');
end

iseries	= sscanf(strseries,'%d');
iexp 	= sscanf(strexp,'%d');
   
%% Do what you have been called to do

switch strarg

    case 'Load' % not really a callback...

        PICK.CurrentlySorting = 'Multispike';
        PICK.vv = {};
        PICK.expt = [];
        [newexpt, nnewrepeats] = ExptLoad(stranimal, iseries, iexp);
        if isempty(newexpt)
            error('Could not load experiment');
        end
        LoadedExpt(newexpt, nnewrepeats, 0); % this SHOULD set PICK.expt

    case 'LoadMichigan' % not really a callback...

        PICK.CurrentlySorting = 'Michigan';
        PICK.vv = {};
        PICK.expt = [];
        [newexpt, nnewrepeats] = ExptLoadMichigan(stranimal, iseries, iexp);
        LoadedExpt(newexpt, nnewrepeats, 0); % this SHOULD set PICK.expt

    case 'LoadCerebusTraces' % not really a callback...

        PICK.CurrentlySorting = 'CerebusTraces';
        PICK.vv = {};
        PICK.expt = [];
        try
            newexpt = ExptLoadCerebusTraces(stranimal, iseries, iexp);
        catch ME
            disp(ME.identifier);
            disp(ME.message);
            disp(ME.cause);
            disp(ME.stack);
            fprintf('Could not load cerebus traces\n');
            return
        end
        if isempty(newexpt)
            fprintf('Could not load cerebus traces\n');
            return
        end
        LoadedExpt(newexpt, newexpt.nrepeats, 0); % this SHOULD set PICK.expt

    case 'LoadContinuousData' % not really a callback...

        PICK.CurrentlySorting = 'CerebusTraces';
        PICK.vv = {};
        PICK.expt = [];
        try
            newexpt = LoadContinuousVRData(stranimal, str2num(strseries), str2num(strexp));
        catch ME
            disp(ME);
            fprintf('Could not load cerebus traces\n');
            return
        end
        if isempty(newexpt)
            fprintf('Could not load cerebus traces\n');
            return
        end
        LoadedExpt(newexpt, 1, 0); % this SHOULD set PICK.expt
        
    case 'pbUpdate'

        switch PICK.CurrentlySorting
            case 'Multispike'
                [newexpt, nnewrepeats] = ExptLoad(PICK.animal, PICK.iseries, PICK.iexp, PICK.expt);
                LoadedExpt(newexpt, nnewrepeats, 1); %
            case 'Michigan'
                [newexpt, nnewrepeats] = ExptLoadMichigan(PICK.animal, PICK.iseries, PICK.iexp, PICK.expt);
                LoadedExpt(newexpt, nnewrepeats, 1); %
            case 'CerebusTraces'
                msgbox('Cannot update cerebus traces -- the experiment is meant to be over');
        end

    case 'chkFilter'

        if get(findobj(PICK.spikesorter,'Tag','chkFilter'),'Value')
            if PICK.expt.samplerate <15000
                % acquired with Multispike or with TDT
                [b,a] = butter(3,2*250/PICK.expt.samplerate,'high');
                fprintf('Highpass Butterworth filter above 250 Hz\n');
            else
                % acquired with traces
                [b,a] = butter(3,2*[250 5000]/PICK.expt.samplerate); % new from 2010-03-16
                fprintf('Bandpass Butterworth filter between 250 Hz and 5 kHz\n');
            end

            % [H,W] = freqz(b,a,512,PICK.expt.samplerate); figure; semilogx(W,abs(H));
            PICK.chans(ichan).filtercoeffs = {b,a};
        else
            PICK.chans(ichan).filtercoeffs = [];
        end

        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        [PICK.chans(ichan), PICK.vv{ichan}] = ChanFindCandidates(PICK.chans(ichan),PICK.expt);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');

        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        PICK.chans(ichan) = ChanFindPrototypes(PICK.vv{ichan},PICK.chans(ichan),irepeat);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');
        if PICK.chans(ichan).nprots==0, fprintf('\n'); end

        FigSpikeSorter_callbacks txtMinCorrelation

    case 'pbChanPlus'

        ichan = ichan+1;
        if ichan>length(PICK.chans) % PICK.expt.nchans
            ichan = 1;
        end
        set(findobj(PICK.spikesorter,'tag','popChannel'),'value',ichan);
        FigSpikeSorter_callbacks popChannel

    case 'pbChanMinus'

        ichan = ichan-1;
        if ichan<1
            ichan = length(PICK.chans); % PICK.expt.nchans;
        end
        set(findobj(PICK.spikesorter,'tag','popChannel'),'value',ichan);
        FigSpikeSorter_callbacks popChannel

    case 'popChannel'

        setfields(ichan);

         % MC 2010-03-26 debugging a strange error
        if MatteoDebug && ~isempty(PICK.chans(ichan).cc) && length(PICK.chans(ichan).spikes)~=length(PICK.chans(ichan).cc)
            error('Weird stuff is happening -- call Matteo');
        end
        
        if ~PICK.chans(ichan).done % IF THIS CHAN IS "DONE" SHOULD NOT RECOMPUTE...

            set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
            [ PICK.chans(ichan), PICK.vv{ichan} ] = ChanFindCandidates(PICK.chans(ichan),PICK.expt);
            set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');

            set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
            PICK.chans(ichan) = ChanFindPrototypes(PICK.vv{ichan},PICK.chans(ichan),irepeat);
            set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');
            if PICK.chans(ichan).nprots==0, fprintf('\n'); end

        end
 
        FigSpikeSorter_callbacks txtMinCorrelation

    case 'txtThreshold'

        thresh = get(findobj(PICK.spikesorter,'Tag','txtThreshold'),'String');
        PICK.chans(ichan).thresh = sscanf(thresh,'%f',1);

         % MC 2010-03-26 debugging a strange error
        if MatteoDebug && ~isempty(PICK.chans(ichan).cc) && length(PICK.chans(ichan).spikes)~=length(PICK.chans(ichan).cc)
            error('Weird stuff is happening -- call Matteo');
        end
        
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        [PICK.chans(ichan), PICK.vv{ichan}] = ChanFindCandidates(PICK.chans(ichan),PICK.expt);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');

        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        PICK.chans(ichan) = ChanFindPrototypes(PICK.vv{ichan},PICK.chans(ichan),irepeat);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');
        if PICK.chans(ichan).nprots==0, fprintf('\n'); end

        % MC 2010-03-26 debugging a strange error
        if MatteoDebug && ~isempty(PICK.chans(ichan).cc) && length(PICK.chans(ichan).spikes)~=length(PICK.chans(ichan).cc)
            error('Weird stuff is happening -- call Matteo');
        end
        
        FigSpikeSorter_callbacks txtMinCorrelation

    case 'txtTStart'

        NewTStart = str2double(get(findobj(PICK.spikesorter,'Tag','txtTStart'),'String'));
        if isnan(NewTStart) || NewTStart <=-3 || NewTStart >= 0
            set(findobj(PICK.spikesorter,'Tag','txtTStart'),'String',-3);
        end

        tlims(1) = str2double(get(findobj(PICK.spikesorter,'Tag','txtTStart'),'String'));
        tlims(2) = str2double(get(findobj(PICK.spikesorter,'Tag','txtTEnd'),'String'));
        set(findobj(PICK.spikesorter,'Tag','axSpikes'),'xlim',tlims);

    case 'txtTEnd'

        NewTEnd = str2double(get(findobj(PICK.spikesorter,'Tag','txtTEnd'),'String'));
        if isnan(NewTEnd) || NewTEnd <=0 || NewTEnd >= 5;
            set(findobj(PICK.spikesorter,'Tag','txtTEnd'),'String',5);
        end

        tlims(1) = str2double(get(findobj(PICK.spikesorter,'Tag','txtTStart'),'String'));
        tlims(2) = str2double(get(findobj(PICK.spikesorter,'Tag','txtTEnd'),'String'));
        set(findobj(PICK.spikesorter,'Tag','axSpikes'),'xlim',tlims);

    case 'txtMinDur'

        mindur = get(findobj(PICK.spikesorter,'Tag','txtMinDur'),'String');
        PICK.chans(ichan).mindur = sscanf(mindur,'%f',1);

        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        [PICK.chans(ichan), PICK.vv{ichan}] = ChanFindCandidates(PICK.chans(ichan),PICK.expt);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');

        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        PICK.chans(ichan) = ChanFindPrototypes(PICK.vv{ichan},PICK.chans(ichan),irepeat);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');
        if PICK.chans(ichan).nprots==0, fprintf('\n'); end

        FigSpikeSorter_callbacks txtMinCorrelation

    case 'rbSignUp'

        if get(findobj(PICK.spikesorter,'tag','rbSignUp'),'Value')==1
            set(findobj(PICK.spikesorter,'tag','rbSignDown'),'Value',0);
            PICK.chans(ichan).threshsign = 1;
        else
            set(findobj(PICK.spikesorter,'tag','rbSignDown'),'Value',1);
            PICK.chans(ichan).threshsign = -1;
        end

        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        [PICK.chans(ichan), PICK.vv{ichan}] = ChanFindCandidates(PICK.chans(ichan),PICK.expt);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');

        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        PICK.chans(ichan) = ChanFindPrototypes(PICK.vv{ichan},PICK.chans(ichan),irepeat);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');
        if PICK.chans(ichan).nprots==0, fprintf('\n'); end

        FigSpikeSorter_callbacks txtMinCorrelation

    case 'rbSignDown'

        if get(findobj(PICK.spikesorter,'tag','rbSignDown'),'Value')==1
            set(findobj(PICK.spikesorter,'tag','rbSignUp'),'Value',0);
            PICK.chans(ichan).threshsign = -1;
        else
            set(findobj(PICK.spikesorter,'tag','rbSignUp'),'Value',1);
            PICK.chans(ichan).threshsign = 1;
        end

        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        [PICK.chans(ichan), PICK.vv{ichan}] = ChanFindCandidates(PICK.chans(ichan),PICK.expt);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');


        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        PICK.chans(ichan) = ChanFindPrototypes(PICK.vv{ichan},PICK.chans(ichan),irepeat);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');
        if PICK.chans(ichan).nprots==0, fprintf('\n'); end

        FigSpikeSorter_callbacks txtMinCorrelation

    case 'popRepeat'

        % if some prototypes are unfrozen (not checked), recompute
        if ~all(PICK.chans(ichan).frozenprots(1:PICK.chans(ichan).nprots))
            PICK.chans(ichan) = ChanFindPrototypes(PICK.vv{ichan},PICK.chans(ichan),irepeat);
        end

        spikeshow(ichan,irepeat);

    case 'txtMinCorrelation'

        % MC 2010-03-26 debugging a strange error
        if MatteoDebug && ~isempty(PICK.chans(ichan).cc) && length(PICK.chans(ichan).spikes)~=length(PICK.chans(ichan).cc)
            error('Weird stuff is happening -- call Matteo');
        end
        
        mincorr = get(findobj(PICK.spikesorter,'Tag','txtMinCorrelation'),'String');
        mincorr = sscanf(mincorr,'%f',1)/100;
        if isempty(mincorr)
            mincorr = 0.5; % changed by MC 2007-09-10
        end
        if mincorr > 1 || mincorr <=0
            set(findobj(PICK.spikesorter,'Tag','txtMinCorrelation'),'String',num2str(100*PICK.chans(ichan).mincorr,3))
            return
        end
        PICK.chans(ichan).mincorr = mincorr;

        %---------------- classify
        if ~isempty(PICK.chans(ichan).spikes)
            PICK.chans(ichan) = ChanSpikeClassify(PICK.vv{ichan},PICK.chans(ichan));
        end

        PICK.chans(ichan).done = 1;
        
        % MC 2010-03-26 debugging a strange error
        if MatteoDebug && ~isempty(PICK.chans(ichan).cc) && length(PICK.chans(ichan).spikes)~=length(PICK.chans(ichan).cc)
            error('Weird stuff is happening -- call Matteo');
        end
        
        if all([PICK.chans.done]) || length(PICK.chans)>16
            set(findobj(PICK.spikesorter,'Tag','pbSave'),'enable','on');
        end

        % the following can be amazingly slow! cla with 15000 traces is super slow...
        spikeshow(ichan, irepeat);

        %-------------- update the controls
        setfields(ichan);

    case 'btnPrototypeON'

        theprot = sscanf(get(gcbo,'UserData'),'%d',1);

        nprots = PICK.chans(ichan).nprots;

        if get(gcbo,'value')
            % ---------- one prototype was added
            % fix nprots
            PICK.chans(ichan).nprots = nprots+1;
        else
            % ---------- one prototype was deleted
            % move up the remaining ones
            for iprot = theprot:(nprots-1)
                PICK.chans(ichan).prots(iprot,:) 		= PICK.chans(ichan).prots(iprot+1,:);
                PICK.chans(ichan).frozenprots(iprot) 	= PICK.chans(ichan).frozenprots(iprot+1);
                PICK.chans(ichan).cellids(iprot) 		= PICK.chans(ichan).cellids(iprot+1);
            end
            % delete the last
            PICK.chans(ichan).prots(nprots,:) = [];
            PICK.chans(ichan).frozenprots(nprots) = 0;
            PICK.chans(ichan).cellids(nprots) = NaN;
            % fix nprots
            PICK.chans(ichan).nprots = nprots-1;
        end

        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
        PICK.chans(ichan) = ChanFindPrototypes(PICK.vv{ichan},PICK.chans(ichan),irepeat);
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');
        if PICK.chans(ichan).nprots==0, fprintf('\n'); end

        FigSpikeSorter_callbacks txtMinCorrelation

    case 'pbSave'

        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');

        if ~all([PICK.chans.done])
            msgbox('You need to look at all channels before you can save');
            success = 0;
        else
            if ~isfield(PICK,'overwriteflag'),
                PICK.overwriteflag = 'ask';
            end
            success = ExptSaveAnalysis(PICK.expt, PICK.chans, DIRS.spikes, DIRS.data, PICK.overwriteflag);
        end

        if success
            set(findobj(PICK.spikesorter,'Tag','pbSave'),'enable','off');
            % send a wake-up call to FigPicker:
            FigPicker_callbacks LoadUnits;
        end
        
        set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');
        
    case 'chkPrototypeFreeze'

        iprot = str2double(get(gcbo,'UserData'));
        PICK.chans(ichan).frozenprots(iprot) = get(gcbo,'value');

        if PICK.chans(ichan).frozenprots(iprot)==0
            % recompute the prototypes
            set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','watch');
            PICK.chans(ichan) = ChanFindPrototypes(PICK.vv{ichan},PICK.chans(ichan),irepeat);
            set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');
            if PICK.chans(ichan).nprots==0, fprintf('\n'); end

            FigSpikeSorter_callbacks txtMinCorrelation

        end

    case 'txtUnitID'

        iprot = str2double(get(gcbo,'UserData'));
        cellid = str2double(get(gcbo,'string'));
        if isempty(cellid)
            set(gcbo,'string',PICK.chans(ichan).cellids(iprot));
            return
        end
        PICK.chans(ichan).cellids(iprot) = cellid;

        set(findobj(PICK.spikesorter,'Tag','pbSave'),'enable','on');

    case 'pbClose'

        celllist = zeros(PICK.chans(ichan).nprots,1);
        for iprot = 1:PICK.chans(ichan).nprots
            celllist(iprot) = PICK.chans(ichan).cellids(iprot); % the cell id associated with the prototype
        end

        set(PICK.spikesorter,'visible','off');

    otherwise

        disp(['Do not know this argument: ' strarg]);
end

set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');


%%

function setfields(ichan)
% SETFIELDS sets all the fields for a given channel
%

global PICK;

set(findobj(PICK.spikesorter,'Tag','txtThreshold'),'String',num2str(PICK.chans(ichan).thresh,4));
set(findobj(PICK.spikesorter,'Tag','txtMinCorrelation'),'String',num2str(100*PICK.chans(ichan).mincorr,3));
if isfield(PICK.chans(ichan),'threshsign') && ~isempty(PICK.chans(ichan).threshsign)
    set(findobj(PICK.spikesorter,'Tag','rbSignUp'  ),'Value',PICK.chans(ichan).threshsign ==  1);
    set(findobj(PICK.spikesorter,'Tag','rbSignDown'),'Value',PICK.chans(ichan).threshsign == -1);
else
    set(findobj(PICK.spikesorter,'Tag','rbSignUp'  ),'Value',1);
    set(findobj(PICK.spikesorter,'Tag','rbSignDown'),'Value',0);
end
if isfield(PICK.chans(ichan),'mindur')
    set(findobj(PICK.spikesorter,'Tag','txtMinDur'  ),'String',PICK.chans(ichan).mindur);
else
    set(findobj(PICK.spikesorter,'Tag','txtMinDur'  ),'String','2');
end

% ---------- prototype info
for iprot = 1:PICK.chans(ichan).nprots
    set(findobj(PICK.spikesorter,'Tag','btnPrototypeON','UserData',num2str(iprot)),...
        'value',1);
    set(findobj(PICK.spikesorter,'Tag','chkPrototypeFreeze','UserData',num2str(iprot)),...
        'value',PICK.chans(ichan).frozenprots(iprot));
    set(findobj(PICK.spikesorter,'Tag','txtUnitID','UserData',num2str(iprot)),...
        'string',num2str(PICK.chans(ichan).cellids(iprot)));
end
for iprot = PICK.chans(ichan).nprots+1:4
    set(findobj(PICK.spikesorter,'Tag','btnPrototypeON','UserData',num2str(iprot)),'value',0);
    set(findobj(PICK.spikesorter,'Tag','txtUnitID','UserData',num2str(iprot)),...
        'string',[]);
    set(findobj(PICK.spikesorter,'Tag','chkPrototypeFreeze','UserData',num2str(iprot)),...
        'value',0);
end

% ----------- filter info
if isfield(PICK.chans(ichan),'filtercoeffs')
    if isempty(PICK.chans(ichan).filtercoeffs)
        set(findobj(PICK.spikesorter,'Tag','chkFilter'),'Value',0);
    else
        set(findobj(PICK.spikesorter,'Tag','chkFilter'),'Value',1);
    end
else
    set(findobj(PICK.spikesorter,'Tag','chkFilter'),'Value',0);
end

%--------- END OF SUBFUNCTION setfields

%%

function plots = spikeshow(ichan, irepeat )
% SPIKESHOW
% plots = spikeshow(ichan, irepeat)
% where plots has fields vv{}, thresh, prots
%
% plots = spikeshow(ichan) shows all repeats

global PICK;

vv = PICK.vv{ichan};
tt = PICK.chans(ichan).tt;

axes(findobj(PICK.spikesorter,'Tag','axSpikes')); cla %#ok<MAXES>
tlims(1) = str2double(get(findobj(PICK.spikesorter,'Tag','txtTStart'),'String'));
tlims(2) = str2double(get(findobj(PICK.spikesorter,'Tag','txtTEnd'),'String'));
set(gca,'xlim',tlims);

if ~isempty(PICK.chans(ichan).spikes)

    if nargin == 1
        % all the repeats
        mask = ones(size(vv,1),1);
    else
        mask = ([PICK.chans(ichan).spikes(:).irpt] == irepeat);
    end
else
    mask = [];
    vv = [];
end

allcolors = 'wrgbcymk';

% title([' Channel ' num2str(ichan)]);

plots.vv = {};
plots.thresh = [];
plots.prots = [];
plots.mindur = [];

hold on

% plot the candidates with the appropriate color
if ~isempty(vv)
    plots.vv = {};
    for iprot = PICK.chans(ichan).nprots:-1:0 % inverted order so unclassified go last
        if any(PICK.chans(ichan).cc==iprot)
            % in the following, you have to use tt(:) or there can be a bug...
            if length(mask) == length(PICK.chans(ichan).cc)
                ppp = find(mask(:) & PICK.chans(ichan).cc==iprot);
            else
                ppp = find( PICK.chans(ichan).cc==iprot);
            end
            if length(ppp)<1
                % do nothing
            elseif length(ppp) > 5000
                % dont' plot all the spikes, it takes forever to delete them!
                fprintf(1,'--------> Plotting a summary of the spike shapes for the %d spikes assigned to prototype %d\n',length(ppp),iprot)
                meanvv = mean(vv(ppp,:),1);
                stdvv = std(vv(ppp,:),0,1);
                plots.vv{iprot+1} = fillplot(tt(:)', meanvv(:)'-3*stdvv(:)', meanvv(:)'+3*stdvv(:)', allcolors(iprot+1));
                mm = repmat(meanvv(:),1,length(ppp))';
                ee = repmat(stdvv(:),1,length(ppp))';
                outliers = find(any(abs(vv(ppp,:)-mm)>3*ee,2));
                if ~isempty(outliers)
                    plots.vv{iprot+1} = [ plots.vv{iprot+1}; plot(tt(:)', vv(outliers,:)', 'color',allcolors(iprot+1))];
                end
                set(plots.vv{iprot+1},'UserData',[ichan, iprot]);

            else
                % you must transpose the vv, otherwise when nrows = ncols...
                plots.vv{iprot+1} = plot(tt(:)', vv(ppp,:)', 'color',allcolors(iprot+1));
                set(plots.vv{iprot+1},'UserData',[ichan, iprot]);
            end
        end
    end
end

% plot the prototypes
ttsamples = find(tt> -0.5 & tt<=PICK.chans(ichan).mindur);
for iprot = 1:PICK.chans(ichan).nprots
    plots.prots(iprot) = plot(tt(ttsamples),PICK.chans(ichan).prots(iprot,ttsamples),'linewidth',2,'color','k');
    set(plots.prots(iprot),'UserData',[ichan, iprot]);
end

plots.thresh = plot([tt(1) tt(end)],PICK.chans(ichan).thresh*[1 1],'-');
set(plots.thresh,'UserData',ichan);

vlim = sort(PICK.chans(ichan).thresh*[-2 2]); % hack

plots.mindur = plot(PICK.chans(ichan).mindur*[1 1], vlim,'--');
set(plots.mindur,'UserData',ichan);

plot([tt(1) tt(end)],vlim,'.'); % just to make sure the range is right

plot([0 0], vlim, ':'); % the y axis at t = 0

set(gca,'ylim',[-inf inf],'xtick',-10:10);
xlabel('Time (ms)');
set(gca,'color',0.5*[1 1 1])

set(cat(1,plots.vv{:}),'buttondownfcn','disp(get(gcbo,''UserData''))');
set(plots.thresh,'buttondownfcn','disp(get(gcbo,''UserData''))');
set(plots.prots,'buttondownfcn','disp([''Prototype '' mat2str(get(gcbo,''UserData''))])');

set(cat(1,plots.vv{:}),'buttondownfcn','BringToTop');
set(plots.thresh,'buttondownfcn','BringToTop');
set(plots.prots,'buttondownfcn','BringToTop');

%--------- END OF SUBFUNCTION spikeshow

%%

function LoadedExpt(newexpt, nnewrepeats, updateflag)
% LoadedExpt -- things to be done after loading or updating an experiment
%

% set(findobj(PICK.spikesorter,'tag','figSpikeSorter'),'pointer','arrow');

if nargin<3, updateflag = 0; end

global PICK DIRS

if isempty(newexpt) || isempty(newexpt.data) || nnewrepeats==0
    return;
end

if ~updateflag
    % loading an exp:
    PICK.overwriteflag = 'ask';
else
    % updating an exp:
    PICK.overwriteflag = 'overwrite';
end

if isfield(PICK,'expt') && isfield(PICK.expt,'animal') && ~strcmp( newexpt.animal, PICK.expt.animal )
    % we have changed animal. SAMPLERATE could be different...
    PICK.chans = [];
end

PICK.expt = newexpt;

for ichan = 1:PICK.expt.nchans % length(PICK.chans)
    PICK.chans(ichan).done = 0;
end

if ~updateflag
    % if loading a new file, initialize the channels
    
    % Matteo moved this code here from old subfunction initialize_chans
    chanloaded = 0;
    
    % see if user wants to retain the discrim parameters
    CanRetain = ~isempty(PICK.chans) && isfield(PICK.chans,'spikes');
    
    fprintf('\nWorking directory: %s\\%s\\%d\n', DIRS.spikes, PICK.expt.animal, PICK.expt.iseries);
    CurrChanFileName = ChanGetFileName(DIRS.spikes,PICK.expt.animal,PICK.expt.iseries,PICK.expt.iexp  ,PICK.CurrentlySorting);
    PrevChanFileName = ChanGetFileName(DIRS.spikes,PICK.expt.animal,PICK.expt.iseries,PICK.expt.iexp-1,PICK.CurrentlySorting);
    NextChanFileName = ChanGetFileName(DIRS.spikes,PICK.expt.animal,PICK.expt.iseries,PICK.expt.iexp+1,PICK.CurrentlySorting);

    CanLoad = false(3,1);
    CanLoad(1) = ~chanloaded && exist(CurrChanFileName,'file')==2; % the current file
    CanLoad(2) = ~chanloaded && exist(PrevChanFileName,'file')==2; % the previous file
    CanLoad(3) = ~chanloaded && exist(NextChanFileName,'file')==2; % the next file
    
    OptionsList = {};
    if CanRetain   , OptionsList{end+1} = 'Retain'; end  
    if any(CanLoad), OptionsList{end+1} = 'Load'; end
    OptionsList{end+1} = 'New';
   
    if numel(OptionsList)>1
        strSource = questdlg('Source of spike discrimination parameters?', 'Spike Sorter', OptionsList{:}, OptionsList{end} );
    else
        strSource = 'New';
    end
     
    if strcmp(strSource, 'Load')
        
        OptionsList = {};
        if CanLoad(1), OptionsList{end+1} = 'Current' ; end
        if CanLoad(2), OptionsList{end+1} = 'Previous'; end
        if CanLoad(3), OptionsList{end+1} = 'Next'    ; end
        
        if numel(OptionsList)>1 || ~CanLoad(1)
            strSource = questdlg('Load discrimination parameters from which experiment?', 'Spike Sorter', OptionsList{:}, OptionsList{1} );
        else
            strSource = 'Current';
        end
    end
    
    switch strSource
        
        case 'Retain'

            PICK.chans = ChanImport(PICK.expt,PICK.chans);
            for ichan = 1:length(PICK.chans)
                PICK.chans(ichan).done = 0;
                PICK.chans(ichan).frozenprots = [1 1 1 1];
                PICK.chans(ichan).spikes = [];
                PICK.chans(ichan).cc = [];
            end

        case 'Current'
            loadedchans = ChanLoad(DIRS.spikes,PICK.expt.animal,PICK.expt.iseries,PICK.expt.iexp,PICK.CurrentlySorting);
            PICK.chans = ChanImport(PICK.expt,loadedchans);
            for ichan = 1:length(PICK.chans)
                PICK.chans(ichan).frozenprots = [1 1 1 1];
                PICK.vv{ichan} = ChanFindCandidates(PICK.chans(ichan),PICK.expt);
            end
        
        case 'Previous'
    
            loadedchans = ChanLoad(DIRS.spikes,PICK.expt.animal,PICK.expt.iseries,PICK.expt.iexp-1,PICK.CurrentlySorting);
            PICK.chans = ChanImport(PICK.expt,loadedchans); % problem here, PICK.expt.channames is 1 2 3 4
            for ichan = 1:length(PICK.chans)
                PICK.chans(ichan).frozenprots = [1 1 1 1];
                PICK.chans(ichan).done = 0;
                PICK.chans(ichan).spikes = [];
                PICK.chans(ichan).cc = [];
            end
            
        case 'Next'
            
            loadedchans = ChanLoad(DIRS.spikes,PICK.expt.animal,PICK.expt.iseries,PICK.expt.iexp+1,PICK.CurrentlySorting);
            PICK.chans = ChanImport(PICK.expt,loadedchans); % problem here, PICK.expt.channames is 1 2 3 4
            for ichan = 1:length(PICK.chans)
                PICK.chans(ichan).frozenprots = [1 1 1 1];
                PICK.chans(ichan).done = 0;
                PICK.chans(ichan).spikes = [];
                PICK.chans(ichan).cc = [];
            end
            
        otherwise
            
            PICK.chans = ChanInitialize(PICK.expt);
            
    end
    
end

set(findobj(PICK.spikesorter,'Tag','txtFilename'),'String',...
    [ PICK.expt.animal ' Exp ' num2str(PICK.expt.iseries) '_' num2str(PICK.expt.iexp) ]);

set(findobj(PICK.spikesorter,'Tag','popChannel'),'String',PICK.expt.channames,'value',1);

if ~updateflag
    % if loading a new file, look at channel one
    set(findobj(PICK.spikesorter,'Tag','popChannel'),'value',1);
end

% enable everything
set(findobj(PICK.spikesorter,'Enable','Off'),'Enable','on');
% disable a couple of controls
if length(PICK.chans) < 2
    set(findobj(PICK.spikesorter,'Tag','popChannel'),'Enable','off');
end

set(findobj(PICK.spikesorter,'Tag','pbSave'),'enable','off');

if ~updateflag
    % if loading a new file, look at repeat one
    set(findobj(PICK.spikesorter,'Tag','popRepeat'),'String',1:PICK.expt.nrepeats,'Value',1)
else
    % go to the last repeat
    set(findobj(PICK.spikesorter,'Tag','popRepeat'),'String',1:PICK.expt.nrepeats,'Value',PICK.expt.nrepeats)
end

FigSpikeSorter_callbacks popChannel

%--------- END OF SUBFUNCTION LoadedExpt


