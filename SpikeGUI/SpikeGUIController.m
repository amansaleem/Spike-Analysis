classdef SpikeGUIController < handle
    
    properties
        figInfo
        DIRS
        animal
        series
        iseries
        CellInfo = cellSelection;
        exptList
        subset
        actionButtons
        data
        numFigObjects = 0;
        FigObjects
        FigObjLisHandles
    end
    
    events 
       UpdateCells
    end
    
    methods
        
        function Controller = SpikeGUIController
            %% Controller figure creation
            SetDirs
            DIRS.spikes = [DIRS.spikes filesep 'Klustered'];
            Controller.DIRS = DIRS;
            createUI(Controller);
        end
        function Controller = createUI(Controller)
            Controller.figInfo.fig = figure('Name', 'Spike Analysis Controller',...
                'Position', [10 300 675 220], ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'NumberTitle', 'off');
            Controller.figInfo.layout = uiextras.HBox('Parent',Controller.figInfo.fig,'Spacing',10,'Padding',5);
            % Column with text and clear button
            Controller.figInfo.textCol = uiextras.VBox('Parent',Controller.figInfo.layout);
            uicontrol('Style','text','Parent',Controller.figInfo.textCol,'String','', 'fontsize',14, 'HorizontalAlignment','left');
            uicontrol('Style','text','Parent',Controller.figInfo.textCol,'String','Animal:', 'fontsize',12, 'HorizontalAlignment','center');
            
            uicontrol('Style','text','Parent',Controller.figInfo.textCol,'String','', 'fontsize',14, 'HorizontalAlignment','left');
            uicontrol('Style','text','Parent',Controller.figInfo.textCol,'String','Series:', 'fontsize',12, 'HorizontalAlignment','center');
            
            uicontrol('Style','text','Parent',Controller.figInfo.textCol,'String','', 'fontsize',14, 'HorizontalAlignment','left');
            uicontrol('Style','text','Parent',Controller.figInfo.textCol,'String','Exp(s):', 'fontsize',12, 'HorizontalAlignment','center');
            
            uicontrol('Style','text','Parent',Controller.figInfo.textCol,'String','', 'fontsize',14, 'HorizontalAlignment','left');
            Controller.figInfo.ClearUI = uicontrol('Parent',Controller.figInfo.textCol,'String','CLEAR', 'fontsize',14, 'Callback', @(src, evt) Controller.clearFunc(Controller));
            
            Controller.figInfo.textCol.Sizes = [25 25 25 25 25 25 25 25] - 3;
            
            % Column with loaders
            Controller.figInfo.loadCol  = uiextras.VBox('Parent',Controller.figInfo.layout,'Spacing',10,'Padding',5);
            Controller.figInfo.animalUI = uicontrol('Parent',Controller.figInfo.loadCol,'String','Select', 'fontsize',12,'Enable','on','Callback',@(~,~) Controller.selectAnimal);
            Controller.figInfo.seriesUI = uicontrol('Parent',Controller.figInfo.loadCol,'String','Select', 'fontsize',12,'Enable','off','Callback',@(~,~) Controller.selectSeries);
            Controller.figInfo.exptUI   = uicontrol('Parent',Controller.figInfo.loadCol,'String','Select', 'fontsize',12,'Enable','off','Callback',@(~,~) Controller.selectExpt);
            
            Controller.figInfo.cellRow   = uiextras.HBox('Parent',Controller.figInfo.loadCol);
            Controller.figInfo.cellMinusUI  = uicontrol('Parent',Controller.figInfo.cellRow,'String','-', 'fontsize',12,'Enable','off','Callback',@(~,~) Controller.subtractCell);
            Controller.figInfo.cellUI       = uicontrol('Style','popup','Parent',Controller.figInfo.cellRow,'String',Controller.CellInfo.CellListString, ...
                'fontsize',12,'Enable','off','Callback',@(src, evt) Controller.selectCell(get(src,'Value')));
            Controller.figInfo.cellPlusUI   = uicontrol('Parent',Controller.figInfo.cellRow,'String','+', 'fontsize',12,'Enable','off','Callback',@(~,~) Controller.addCell);
            Controller.figInfo.cellRow.Sizes = [-1 -3 -1];
            set(Controller.figInfo.cellRow,'Spacing',5);
            
            Controller.figInfo.loadCol.Sizes = [40 40 40 40];
            
            
            % Column with subset choice
%             Controller.figInfo.subsetCol = uiextras.VBox('Parent',Controller.figInfo.layout,'Spacing',10);
%             Controller.figInfo.subsetCol2 = uiextras.Grid('Parent',Controller.figInfo.subsetCol,'Spacing',10);
%             Controller.figInfo.ResetUI = uicontrol('Parent',Controller.figInfo.subsetCol, 'String','RESET', 'fontsize',14);
%             Controller.figInfo.subsetCol.Sizes = [135 45];
%             
%             uicontrol('Style','text','Parent',Controller.figInfo.subsetCol2,'String','ROOM', 'fontsize',12);
%             Controller.figInfo.ContrastUI = uicontrol('Style','popup','Parent',Controller.figInfo.subsetCol2,'String','Contrast', 'fontsize',12);
%             Controller.figInfo.GainUI = uicontrol('Style','popup','Parent',Controller.figInfo.subsetCol2,'String','Gain', 'fontsize',12);
%             Controller.figInfo.LengthUI = uicontrol('Style','popup','Parent',Controller.figInfo.subsetCol2,'String','RoomLength', 'fontsize',12);
%             
%             uicontrol('Style','text','Parent',Controller.figInfo.subsetCol2,'String','REWARD', 'fontsize',12);
%             Controller.figInfo.StateUI = uicontrol('Style','popup','Parent',Controller.figInfo.subsetCol2,'String','State', 'fontsize',12);
%             Controller.figInfo.ResultUI = uicontrol('Style','popup','Parent',Controller.figInfo.subsetCol2,'String','Result', 'fontsize',12);
%             Controller.figInfo.LocationUI = uicontrol('Style','popup','Parent',Controller.figInfo.subsetCol2,'String','Location', 'fontsize',12);
%             
%             set(Controller.figInfo.subsetCol2,'ColumnSizes',[-1 -1],'RowSizes',[-1 -1 -1 -1])
            
            %Column with run boxes
            Controller.figInfo.buttonCol = uiextras.VBox('Parent',Controller.figInfo.layout,'Spacing',10, 'Padding',5);
            Controller.figInfo.BehSpksUI = uicontrol('Parent',Controller.figInfo.buttonCol,'String','Behav. + Spks', 'fontsize',14,'callback',@(~,~) Controller.startBehSpk, 'Enable','off');
            Controller.figInfo.Map1DUI = uicontrol('Parent',Controller.figInfo.buttonCol,'String','Map 1D', 'fontsize',14,'callback',@(~,~) Controller.startMap1D, 'Enable','off');
            Controller.figInfo.Map2DUI = uicontrol('Parent',Controller.figInfo.buttonCol,'String','Map 2D', 'fontsize',14,'callback',@(~,~) Controller.startMap2D, 'Enable','off');
            
            Controller.figInfo.buttonCol.Sizes = [50 50 50];
            
            % Set sizes
            Controller.figInfo.layout.Sizes = [-1 -2 -1];
        end
        %% Callbacks
        function startBehSpk(obj)
            obj.numFigObjects = obj.numFigObjects + 1;
            eval(['obj.FigObjects.figGUIobj' num2str(obj.numFigObjects) ' = BehavSpkGUI(obj.data.es);']);

            obj.FigObjLisHandles{obj.numFigObjects} = obj.addlistener('UpdateCells',eval(['@(~,~) obj.FigObjects.figGUIobj' num2str(obj.numFigObjects) '.makeAllPlots(obj.CellInfo.ChosenCell);']));
            eval(['obj.FigObjects.figGUIobj' num2str(obj.numFigObjects) '.makeAllPlots(obj.CellInfo.ChosenCell);']);
        end
        
        function startMap1D(obj)
            obj.numFigObjects = obj.numFigObjects + 1;
            eval(['obj.FigObjects.figGUIobj' num2str(obj.numFigObjects) ' = Map1DGUI(obj.data.es);']);

            obj.FigObjLisHandles{obj.numFigObjects} = obj.addlistener('UpdateCells',eval(['@(~,~) obj.FigObjects.figGUIobj' num2str(obj.numFigObjects) '.makeAllPlots(obj.CellInfo.ChosenCell);']));
            eval(['obj.FigObjects.figGUIobj' num2str(obj.numFigObjects) '.makeAllPlots(obj.CellInfo.ChosenCell);']);
        end
        function startMap2D(obj)
            obj.numFigObjects = obj.numFigObjects + 1;
            eval(['obj.FigObjects.figGUIobj' num2str(obj.numFigObjects) ' = Map2DGUI(obj.data.es);']);

            obj.FigObjLisHandles{obj.numFigObjects} = obj.addlistener('UpdateCells',eval(['@(~,~) obj.FigObjects.figGUIobj' num2str(obj.numFigObjects) '.makeAllPlots(obj.CellInfo.ChosenCell);']));
            eval(['obj.FigObjects.figGUIobj' num2str(obj.numFigObjects) '.makeAllPlots(obj.CellInfo.ChosenCell);']);
        end
        function selectAnimal(obj)
            dir_list = dir(obj.DIRS.ball);
            animal_list = [];
            for n = 3:length(dir_list)
                mdate(n-2) = dir_list(n).datenum;
            end
            [~, dorder] = sort(mdate);
            for n = length(dorder):-1:1
                pos = strfind(dir_list(dorder(n)+2).name,'_BALL');
                if length(pos)>0
                    animal_list = [animal_list, {dir_list(dorder(n)+2).name}];
                end
                %%%% Adding the next three lines
                pos = strfind(dir_list(dorder(n)+2).name,'JF');
                if length(pos)>0
                    animal_list = [animal_list, {dir_list(dorder(n)+2).name}];
                end
            end
            
            [Selection,ok] = listdlg('PromptString','Select an animal:',...
                'SelectionMode','single',...
                'ListString',animal_list);
            if ok
                obj.animal = animal_list(Selection);
                obj.animal = obj.animal{1};
                set(obj.figInfo.seriesUI, 'Enable', 'on');
                set(obj.figInfo.animalUI, 'Enable', 'off','String',obj.animal);
            end
        end
        function selectSeries(obj)
            infoAll = getDataInfo(obj.animal);
            date_list = [];
            for n = 1:length(infoAll)
                date_list = [date_list, {infoAll(n).date}];
            end
            
            if ismac
                date_list{1} = 'bloddy mac';
            end
            [Selection,ok] = listdlg('PromptString','Select a series:',...
                'SelectionMode','single',...
                'ListString',date_list);
            if ok
                date = date_list(Selection);
                obj.series = date{1};
                obj.iseries = str2num(date{1});
                
                set(obj.figInfo.exptUI, 'Enable', 'on');
                set(obj.figInfo.seriesUI, 'Enable', 'off','String',obj.series);
            end
        end     
        function selectExpt(obj)
            infoAll = getDataInfo(obj.animal);
            date_list = [];
            for n = 1:length(infoAll)
                date_list = [date_list, {infoAll(n).date}];
            end
            
            whichdate = obj.series;
            for i=1:numel(infoAll)
                if strcmpi(infoAll(i).date,whichdate)==1
                    expt_info = infoAll(i);
                    break;
                end
            end
            
            for iexp = 1:length(expt_info.sessions)
                expt_list{iexp} = num2str(expt_info.sessions(iexp));
            end
            
            [Selection,ok] = listdlg('PromptString','Select experiment(s):',...
                'SelectionMode','multiple',...
                'ListString',expt_list);
            if ok
                obj.exptList = expt_info.sessions(Selection);
                set(obj.figInfo.exptUI, 'Enable', 'off','String',num2str(obj.exptList));
            end
            try
                obj.data.es = VRLoadMultipleExpts(obj.animal, obj.iseries, obj.exptList,'SPIKES');
                obj.data.es = getESDataSubset(obj.data.es);
                set(obj.figInfo.cellMinusUI,'Enable','on');
                obj.CellInfo.defineCellInfo(obj.data.es.spikeIDs);
                
                set(obj.figInfo.cellUI     ,'Enable','on', 'String',obj.CellInfo.CellListString);
                set(obj.figInfo.cellPlusUI ,'Enable','on');
                set(obj.figInfo.BehSpksUI, 'Enable','on');
                set(obj.figInfo.Map1DUI, 'Enable','on');
                set(obj.figInfo.Map2DUI, 'Enable','on');
            catch
                errordlg('No Spiking data available', 'WARNING!');
                obj.data.es = VRLoadMultipleExpts(obj.animal, obj.iseries, obj.exptList);
            end
        end
        function clearFunc(obj,~,~)
            obj.animal = [];
            obj.series = [];
            obj.exptList = [];
            
            set(obj.figInfo.animalUI, 'Enable', 'on','String','Select');
            set(obj.figInfo.exptUI, 'Enable', 'off','String','Select');
            set(obj.figInfo.seriesUI, 'Enable', 'off','String','Select');
            
            for ifig = 1:obj.numFigObjects
                delete(obj.FigObjLisHandles{ifig});
                name = ['figGUIobj' num2str(ifig)];
                eval(['close(obj.FigObjects.figGUIobj' num2str(ifig) '.FigInfo.fig)']);
                rmfield(obj.FigObjects, name);
            end
            obj.numFigObjects = 0;
        end
        function addCell(obj)
            cellPlus(obj.CellInfo);
            set(obj.figInfo.cellUI, 'Value', obj.CellInfo.ChosenCell)
            notify(obj,'UpdateCells'); 
        end
        function subtractCell(obj)
            cellMinus(obj.CellInfo);
            set(obj.figInfo.cellUI, 'Value', obj.CellInfo.ChosenCell)
            notify(obj,'UpdateCells'); 
        end
        function obj = selectCell(obj, val)
            set(obj.CellInfo, 'ChosenCell', val);
            notify(obj,'UpdateCells'); 
        end
    end
end