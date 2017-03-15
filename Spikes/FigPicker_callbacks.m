function FigPicker_callbacks (strarg)
% FigPicker_callbacks is a switchyard for the FigPicker figure
% 2009-01 AZ modified for new DirDialog function, which gives more options
%    when choosing data directory.  Allows user to cancel dir selection.
% 2009-02 AZ now DirDialog function will be skipped if spikes.m passes
%    along which specific animal the user requested figpicker to open (in
%    PICK)
% 2009-04 MC modified case when there is a single series
% 2009-05 MC added choice of unit
% 2010-05 MC introduced ProtocolDescribe
% 2010-06 MC avoid crash with experiments with no protocol 
% 2010-06 MC made big changes to allow analysis of traces

global PICK DIRS;
persistent units;

pickerfig = PICK.picker; % used to be gcbf...

if nargin<1
    error('FigPicker_callbacks must have one argument');
end

if ~isfield(DIRS,'data'), 		DIRS.data 	= ''; end
if ~isfield(DIRS,'spikes'), 	DIRS.spikes = ''; end

if ~isfield(PICK,'animal'), 	PICK.animal 	= ''; end
if ~isfield(PICK,'iseries'), 	PICK.iseries 	= []; end
if ~isfield(PICK,'iexp'), 		PICK.iexp 		= []; end

set(pickerfig,'name','ExptPick -- Experiment Picker');
set(pickerfig,'pointer','watch'); drawnow

%% added by Matteo 2010-12-13: check that FigSpikeSorter still exists

if ~isfield(PICK,'spikesorter') || ~any(PICK.spikesorter==get(0,'children'))
	PICK.spikesorter = FigSpikeSorter;
    set(PICK.spikesorter,...
        'HandleVisibility','callback', ...
        'Visible','Off',...
        'WindowStyle','Normal',...
        'RendererMode','auto',...
        'resize','off'); % 'CloseRequestFcn',@DenyInteractiveCloseRequest);
end

%% Main loop

switch strarg
    
    case 'pbAnimal'
        
        % olddatadir = DIRS.data;
        
        % this code should call ListAnimals
        if ~isfield(PICK,'exptinfos') && ~isempty(PICK.animal)
            % If an animal is already 'picked'/specified by spikes.m, set it, and
            % skip datadir dialog box
            datadir = DIRS.data;
            strarrAnimals{1} = PICK.animal;
        else
            [datadir,strarrAnimals] = DirDialog('data');
        end
        
        if datadir
            if isdir(datadir)
                
                if ~isempty(strarrAnimals)
                    DIRS.data = datadir;
                    
                    strarrAnimals = sort(upper(strarrAnimals));
                    set(findobj(pickerfig,'tag','popAnimal'),'String',strarrAnimals,'Value',length(strarrAnimals),'enable','on');
                    
                    if ~isfield(PICK,'exptinfos') && ~isempty(PICK.animal) && ~isempty(PICK.iseries)
                        if exist(fullfile(DIRS.data,PICK.animal),'dir')~=7
                            set(pickerfig,'pointer','arrow');
                            errordlg(['Directory ' fullfile(DIRS.data,PICK.animal) ' does not exist or is not readable'], 'Error');
                            return
                        end
                        
                        set(findobj(pickerfig,'enable','on'),'enable','off');
                        set(findobj(pickerfig,'Tag','pbAnimal'),'enable','on');
                        set(findobj(pickerfig,'Tag','popAnimal'),'enable','on');
                        set(findobj(pickerfig,'Tag','popSeries'),'String',num2str(PICK.iseries),'value',1,'enable','on');
                        
                        PICK.exptinfos = ExptReadLogFile(DIRS.data,PICK.animal);
                        
                        if ~isempty(PICK.iexp)
                            
                            protocol = ProtocolLoad(PICK.animal,PICK.iseries,PICK.iexp);
                            strexplist{1} = [ num2str(PICK.iexp), ' -- ', protocol.description ];
                            
                            set(findobj(pickerfig,'String','Series'),'enable','on'); % this is just a label, no reason to disable it.
                            set(findobj(pickerfig,'Tag','popExperiment'),'String',strexplist,'value',1,'enable','on');
                            
                            FigPicker_callbacks popExperiment;
                            set(pickerfig,'pointer','arrow');
                            return;
                        end
                        
                        FigPicker_callbacks popSeries;
                        set(pickerfig,'pointer','arrow');
                        return;
                    end
                    
                    FigPicker_callbacks popAnimal;
                    set(pickerfig,'pointer','arrow');
                    return;
                    
                else
                    warning('Could not find any animals in this directory.  Please choose another directory.')
                    set(pickerfig,'pointer','arrow');
                    return;
                end
            else
                warning([datadir,' is not a valid directory.  Please choose another directory.'])
                set(pickerfig,'pointer','arrow');
                return;
            end
        else
            disp('No directory selected.  Please try again.')
            set(pickerfig,'pointer','arrow');
            return;
        end
        
    case 'popAnimal'
        
        strarrAnimals    = get(findobj(pickerfig,'tag','popAnimal'),'String');
        ianimal          = get(findobj(pickerfig,'tag','popAnimal'),'Value');
        animal           = strarrAnimals{ianimal};
        
        if exist(fullfile(DIRS.data,animal),'dir')~=7
            set(pickerfig,'pointer','arrow');
            errordlg(['Directory ' fullfile(DIRS.data,animal) ' does not exist or is not readable'], 'Error');
            return
        end
        
        serieslist = ListSeries(DIRS.data, animal);
        if isempty(serieslist),
            set(findobj(pickerfig,'enable','on'),'enable','off');
            set(findobj(pickerfig,'Tag','pbAnimal'),'enable','on');
            set(findobj(pickerfig,'Tag','popAnimal'),'enable','on');
           
            set(findobj(pickerfig,'Tag','popSeries'),'String',{''},'value',1);
            set(findobj(pickerfig,'Tag','popExperiment'),'String',{''},'value',1);
            set(findobj(pickerfig,'Tag','listParameters'),'String',{''},'value',1);
            set(findobj(pickerfig,'Tag','listComments'),'String',{'No data'},'value',1);
            
            set(pickerfig,'pointer','arrow');
            return;
        end
        
        % convert to array of strings
        strarrSeries = cell(length(serieslist),1);
        for is = 1:length(serieslist)
            strarrSeries{is} = num2str(serieslist(is));
        end
        
        set(findobj(pickerfig,'enable','on'),'enable','off');
        set(findobj(pickerfig,'Tag','pbAnimal'),'enable','on');
        set(findobj(pickerfig,'Tag','popAnimal'),'enable','on');
        
        set(findobj(pickerfig,'Tag','popSeries'),'String',strarrSeries,'value',1,'enable','on');
        
        PICK.animal = animal;
        PICK.exptinfos = ExptReadLogFile(DIRS.data,animal);
        
        FigPicker_callbacks popSeries;
        
    case 'popSeries'
        
        serieslist    =            get(findobj(pickerfig,'Tag','popSeries'),'String') ;
        strseries     = serieslist{get(findobj(pickerfig,'Tag','popSeries'),'Value' )};
        
        Exptlist = ListExpts(DIRS.data, PICK.animal, str2double(strseries));
        
        % convert to an array of strings
        strarrExpt = cell(length(Exptlist),1);
        for ie = 1:length(Exptlist)
            strarrExpt{ie} = num2str(Exptlist(ie));
        end
        
        PICK.iseries = sscanf(strseries,'%d');
        
        if isempty(strarrExpt),
            set(pickerfig,'pointer','arrow');
            return;
        end
        
        strexplist = cell(length(strarrExpt),1);
        for ie = 1:length(strarrExpt)
            iexp = sscanf(strarrExpt{ie},'%d');
            protocol = ProtocolLoad(PICK.animal,PICK.iseries,iexp);
            if isfield(protocol,'description')
                strexplist{ie} = [ strarrExpt{ie}, ' -- ', protocol.description ];
            else
                strexplist{ie} = [ strarrExpt{ie}, ' -- Missing protocol data' ];
            end
        end
        
        set(findobj(pickerfig,'enable','on'),'enable','off');
        set(findobj(pickerfig,'String','Series'),'enable','on'); % this is just a label, no reason to disable it.
        
        % AZ 2009-02-26: disable pushbutton if no new series in animal
        if str2double(strseries) < length(serieslist)
            set(findobj(pickerfig,'Tag','pbNextSeries'),'enable','on');
        end
        set(findobj(pickerfig,'Tag','pbAnimal'),'enable','on');
        set(findobj(pickerfig,'Tag','popAnimal'),'enable','on');
        set(findobj(pickerfig,'Tag','popSeries'),'enable','on');
        set(findobj(pickerfig,'Tag','popExperiment'),'String',strexplist,'value',1,'enable','on');
        
        FigPicker_callbacks popExperiment;
        
    case 'pbNextSeries'
        
        serieslist = get(findobj(pickerfig,'Tag','popSeries'),'String');
        nseries = length(serieslist);
        
        iserieslist = get(findobj(pickerfig,'Tag','popSeries'),'Value');
        
        if iserieslist<nseries
            % show the next series
            set(findobj(pickerfig,'Tag','popSeries'),'Value',iserieslist+1)
            FigPicker_callbacks popSeries;
        else
            % see if you can load more experiments
            FigPicker_callbacks popAnimal
            if length(get(findobj(pickerfig,'Tag','popSeries'),'String'))>nseries
                set(findobj(pickerfig,'Tag','popSeries'),'Value',iserieslist+1)
                FigPicker_callbacks popSeries;
            else
                set(findobj(pickerfig,'Tag','popSeries'),'Value',iserieslist)
                FigPicker_callbacks popSeries;
            end
        end
        
    case 'popExperiment'
        
        iexplist = get(findobj(pickerfig,'Tag','popExperiment'),'Value');
        explist = get(findobj(pickerfig,'Tag','popExperiment'),'String');
        
        iexp = sscanf(explist{iexplist},'%d ',1);
        
        PICK.iexp = iexp;
        PICK.protocol = [];  % this will be filled up shortly
        PICK.expt = [];      % this will be filled up if user click on "Spikes"
        
        % AZ 2009-02-26: disable pushbutton if no new exp in series
        % MC 2010-05-04: commented out: a new exp may have been run in the meantime
%         if iexplist < length(explist)
             set(findobj(pickerfig,'Tag','pbNextExperiment'),'enable','on');
%         else
%             set(findobj(pickerfig,'Tag','pbNextExperiment'),'enable','off');
%         end
        
        %% DESCRIBE THE PARAMETERS
        protocol = ProtocolLoad( PICK.animal, PICK.iseries, PICK.iexp );
        
        ppp = ProtocolDescribe(protocol);
        set(findobj(pickerfig,'Tag','listParameters'),'String',ppp,'value',1,'enable','on');
        
        %% enable the buttons that inspect parameters, gratings, and show movies
        set(findobj(pickerfig,'Tag','pbInspect'),'enable','on');
        set(findobj(pickerfig,'Tag','pbGratings'),'enable','on');
        set(findobj(pickerfig,'Tag','pbWatch'),'enable','on');
        
        PICK.protocol = protocol;
        
        %% DISPLAY THE COMMENTS
        
        comment = {'',''};
        iexptinfo = find([PICK.exptinfos(:).iseries]==PICK.iseries & [PICK.exptinfos(:).iexp]==PICK.iexp );
        if length(iexptinfo)==1
            exptinfo = PICK.exptinfos(iexptinfo);
            % we get a crash here is the format of the date is
            % DD/MM/YYYY...
            % works best when the format is: '26-Jul-00 13:52'
            comment{1} = [ datestr(datenum(exptinfo.StartTime),'HH:MM') ' ' exptinfo.StartComment ];
            if isempty(exptinfo.EndTime) % for example, if the expt is in progress
                comment(2) = '';
            else
                comment{2} = [ datestr(datenum(exptinfo.EndTime  ),'HH:MM') ' ' exptinfo.EndComment ];
            end
        end
        set(findobj(pickerfig,'Tag','listComments'),'String',comment,'enable','on');
        
        %% set up the appropriate buttons for spike sorting
        
        DataTypes = ProtocolGetDataTypes(PICK);
        p = PICK.protocol;   % I assume we don't need another ProtocolLoad    
        
        % start pessimistic
        set(findobj(pickerfig,'Tag','pbMultispike')       ,'Enable','off');
        set(findobj(pickerfig,'Tag','pbTDT')              ,'Enable','off');
        set(findobj(pickerfig,'Tag','pbCerebusSnippets')  ,'Enable','off');
        set(findobj(pickerfig,'Tag','pbCerebusTraces')    ,'Enable','off');
        set(findobj(pickerfig,'Tag','pbTraceInspectRaw')     ,'Enable','off');
                
        if ~isempty(p) && isfield(p,'nrepeats') && p.nrepeats>0
            if  DataTypes.Traces
                set(findobj(pickerfig,'Tag','pbMultispike')  ,'Enable','on');
                set(findobj(pickerfig,'Tag','pbTraceInspectRaw'),'Enable','on');
            end
            
            if DataTypes.TDT
                set(findobj(pickerfig,'Tag','pbTDT'),'Enable','on');
            end
            
            if DataTypes.CerebusSnippets
                set(findobj(pickerfig,'Tag','pbCerebusSnippets'),'Enable','on');
            end
            
            if DataTypes.CerebusTraces
                set(findobj(pickerfig,'Tag','pbCerebusTraces'),'Enable','on');
                set(findobj(pickerfig,'Tag','pbTraceInspectRaw'),'Enable','on');
            end
        end
        
        FigPicker_callbacks LoadUnits
        
    case 'LoadUnits'
        
        units = [];
        if isfield(DIRS,'spikes')
            % this is a persistent variable
            units = UnitLoad( DIRS.spikes, PICK.animal, PICK.iseries, PICK.iexp, 'batch' );
        end
        
        % the controls that let you work on spike data
        ExptControls = [
            findobj(pickerfig,'Tag','popUnit');
            findobj(pickerfig,'Tag','pbTuning');
            findobj(pickerfig,'Tag','pbRasters');
            findobj(pickerfig,'Tag','pbHistos');
            findobj(pickerfig,'Tag','pbFreq');
            findobj(pickerfig,'Tag','pbCycles')];
        
        popUnit = findobj(pickerfig,'Tag','popUnit');
        set(popUnit,'Value',1);
        
        if length(units)>=1
            set(ExptControls,'Enable','on');
            UnitList = cell(length(units)+1,1);
            UnitList{1} = 'All';
            for iu = 1:length(units)
                UnitList{iu+1} = sprintf('%3d: %s',iu,units(iu).id );
            end
            set(popUnit,'String',UnitList);
            
            % set(popUnit,'string',union({'All'},{units.id}));
        else
            set(ExptControls,'Enable','off');
            set(popUnit,'string','None');
        end
        
    case 'pbNextExperiment'
        
        explist = get(findobj(pickerfig,'Tag','popExperiment'),'String');
        nexps = length(explist);
        
        iexplist = get(findobj(pickerfig,'Tag','popExperiment'),'Value');
        
        if iexplist<nexps
            % show the next experiment
            set(findobj(pickerfig,'Tag','popExperiment'),'Value',iexplist+1)
            FigPicker_callbacks popExperiment;
        else
            % see if you can load more experiments
            FigPicker_callbacks popSeries
            if length(get(findobj(pickerfig,'Tag','popExperiment'),'String'))>nexps
                set(findobj(pickerfig,'Tag','popExperiment'),'Value',iexplist+1)
                FigPicker_callbacks popExperiment;
            else
                set(findobj(pickerfig,'Tag','popExperiment'),'Value',iexplist)
                FigPicker_callbacks popExperiment;
            end
        end
        
    case 'pbMultispike'

        % set the values of animal, series, and experiment on figSpikeSorter
        set( findobj(PICK.spikesorter,'Tag','txtAnimal'), 		'String', PICK.animal );
        set( findobj(PICK.spikesorter,'Tag','txtSeries'), 		'String', num2str(PICK.iseries) );
        set( findobj(PICK.spikesorter,'Tag','txtExperiment'),	'String', num2str(PICK.iexp) );
        
        set(PICK.spikesorter,'Visible','On');
        
        % set(PICK.spikesorter,'HandleVisibility','On');
        FigSpikeSorter_callbacks Load
        % set(PICK.spikesorter,'HandleVisibility','Callback');
        
    case 'pbTDT'
        
        % set the values of animal, series, and experiment on figSpikeSorter
        set( findobj(PICK.spikesorter,'Tag','txtAnimal'), 		'String', PICK.animal );
        set( findobj(PICK.spikesorter,'Tag','txtSeries'), 		'String', num2str(PICK.iseries) );
        set( findobj(PICK.spikesorter,'Tag','txtExperiment'),	'String', num2str(PICK.iexp) );
        
        set(PICK.spikesorter,'Visible','On');
        
        % set(PICK.spikesorter,'HandleVisibility','On');
        FigSpikeSorter_callbacks LoadMichigan
        % set(PICK.spikesorter,'HandleVisibility','Callback');
        
    case 'pbCerebusSnippets'
        
      %    sortnev2(PICK.animal, PICK.iseries, PICK.iexp);
      % explorenev6(PICK.animal, PICK.iseries, PICK.iexp);
        explorenev8(PICK.animal, PICK.iseries, PICK.iexp);
        
    case 'pbCerebusTraces'
        
        % set the values of animal, series, and experiment on figSpikeSorter
        set( findobj(PICK.spikesorter,'Tag','txtAnimal'), 		'String', PICK.animal );
        set( findobj(PICK.spikesorter,'Tag','txtSeries'), 		'String', num2str(PICK.iseries) );
        set( findobj(PICK.spikesorter,'Tag','txtExperiment'),	'String', num2str(PICK.iexp) );
        
        set(PICK.spikesorter,'Visible','On');
        
        % set(PICK.spikesorter,'HandleVisibility','On');
        FigSpikeSorter_callbacks LoadCerebusTraces
        % set(PICK.spikesorter,'HandleVisibility','Callback');

    case 'pbTraceInspectRaw'
             
        TraceInspector(PICK.protocol);
        
    case 'pbTuning'
        
        protocol = ProtocolLoad( PICK.animal, PICK.iseries, PICK.iexp, 'loadscreen' );
        
        graphinfo = ProtocolGetGraphInfo(protocol);
        
        if isempty(graphinfo) % e.g. if the user hits Cancel
            set(pickerfig,'pointer','arrow');
            return
        end
        
        % no need to load units, it is persistent
        % units = UnitLoad( DIRS.spikes, PICK.animal, PICK.iseries, PICK.iexp );
        
        AnalysisType=questdlg('What response measure do you want to plot?', ...
            'Spike Sorter', ...
            'Mean','First Harmonic', 'Second Harmonic', 'Mean');
        
        popUnit = findobj(pickerfig,'Tag','popUnit');
        if get(popUnit,'Value') == 1
            % show results for all units
            for iunit = 1:length(units)
                UnitPlotTuning( units(iunit), protocol, graphinfo, [], AnalysisType );
            end
        else
            % show results only for one unit
            iunit = get(popUnit,'Value')-1;
            UnitPlotTuning( units(iunit), protocol, graphinfo, [], AnalysisType );
        end
        
    case 'pbRasters'
        
        % no need to load units, it is persistent
        % units = UnitLoad( DIRS.spikes, PICK.animal, PICK.iseries, PICK.iexp );
        protocol = ProtocolLoad( PICK.animal, PICK.iseries, PICK.iexp );
        
        popUnit = findobj(pickerfig,'Tag','popUnit');
        if get(popUnit,'Value') == 1
            % show results for all units
            for iunit = 1:length(units)
                UnitPlotRasters( units(iunit), protocol );
            end
        else
            % show results only for one unit
            iunit = get(popUnit,'Value')-1;
            UnitPlotRasters( units(iunit), protocol );
        end
        
    case 'pbHistos'
        
        popUnit = findobj(pickerfig,'Tag','popUnit');
        if get(popUnit,'Value') == 1
            % show results for all units
            for iunit = 1:length(units)
                UnitPlotHistos( units(iunit) );
            end
        else
            % show results only for one unit
            iunit = get(popUnit,'Value')-1;
            UnitPlotHistos( units(iunit) );
        end
        
    case 'pbFreq'
        
        popUnit = findobj(pickerfig,'Tag','popUnit');
        if get(popUnit,'Value') == 1
            % show results for all units
            for iunit = 1:length(units)
                UnitPlotFreq( units(iunit) );
            end
        else
            % show results only for one unit
            iunit = get(popUnit,'Value')-1;
            UnitPlotFreq( units(iunit) );
        end
        
    case 'pbCycles'
        
        
        % to decide the scale
        protocol = ProtocolLoad( PICK.animal, PICK.iseries, PICK.iexp);
        minperiod = min(1./protocol.pfilefreqs);
        
        popUnit = findobj(pickerfig,'Tag','popUnit');
        if get(popUnit,'Value') == 1
            % show results for all units
            for iunit = 1:length(units)
                UnitPlotCycles( units(iunit), minperiod/20 );
            end
        else
            % show results only for one unit
            iunit = get(popUnit,'Value')-1;
            UnitPlotCycles( units(iunit), minperiod/20 );
        end
        
    case 'listParameters'
        
        % do nothing
        
    otherwise
        
        disp(['Do not know argument ' strarg]);
end

set(pickerfig,'pointer','arrow');
