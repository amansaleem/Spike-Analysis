classdef BehavOnlyGUI < SpikeDisplayGUI
    
    properties
        variable
    end
    
    methods
        function obj = BehavOnlyGUI(data)
            obj.Data.es = data;
            obj = obj.CalculateSubsets;
            
            obj.FigInfo.fig = figure('Name', 'Behaviour Only Display',...
                'Position', [50 300 1000 800], ...
                'NumberTitle', 'off');
            obj.FigInfo.layout = uiextras.VBox('Parent', obj.FigInfo.fig);
            obj.FigInfo.titleBar = uiextras.HBox('Parent', obj.FigInfo.layout);
            obj.FigInfo.FigGrid = uiextras.Grid('Parent', obj.FigInfo.layout);
            obj.FigInfo.layout.Sizes = [80 -1];
            set(obj.FigInfo.titleBar, 'Spacing',20)
            
            obj.FigInfo.cellSection = uiextras.VBox('Parent', obj.FigInfo.titleBar);
            obj.FigInfo.cellDisplay = uicontrol('Style','text','Parent', obj.FigInfo.cellSection,'fontsize',12,...
                'String',['Cell ' num2str(obj.ChosenCell) ': ' obj.Data.es.spikeIDs{obj.ChosenCell}]);
            
            obj.FigInfo.cellSelection  = uiextras.HBox('Parent', obj.FigInfo.cellSection);
            obj.FigInfo.cellMinus   = uicontrol('String','-','Parent', obj.FigInfo.cellSelection,'Callback',@obj.cellMinus);
            obj.FigInfo.cellChoose  = uicontrol('Style','popup','String','Cell','Parent', obj.FigInfo.cellSelection,'Callback',@obj.cellChoose);
            obj.FigInfo.cellPlus    = uicontrol('String','+','Parent', obj.FigInfo.cellSelection,'Callback',@obj.cellPlus);
            obj.FigInfo.cellSelection.Sizes = [-1 -3 -1];
            
            obj.FigInfo.cellSection.Sizes = [-1 -1];
            
            obj.FigInfo.varChooseGrid = uiextras.Grid('Parent',obj.FigInfo.titleBar);
            obj.FigInfo.varTraj = uicontrol('Enable','off','Parent',obj.FigInfo.varChooseGrid,'String','Traj','Callback',@obj.trajChoose);
            obj.FigInfo.varVisSpd = uicontrol('Enable','on','Parent',obj.FigInfo.varChooseGrid,'String','Vis Spd','Callback',@obj.visSpdChoose);
            obj.FigInfo.varRunSpd = uicontrol('Enable','on','Parent',obj.FigInfo.varChooseGrid,'String','Run Spd','Callback',@obj.runSpdChoose);
            obj.FigInfo.varDist = uicontrol('Enable','on','Parent',obj.FigInfo.varChooseGrid,'String','Distance','Callback',@obj.distChoose);
            obj.FigInfo.varTotDist = uicontrol('Enable','on','Parent',obj.FigInfo.varChooseGrid,'String','Total Distance','Callback',@obj.totDistChoose);
            obj.FigInfo.varRoomPercent = uicontrol('Enable','on','Parent',obj.FigInfo.varChooseGrid,'String','Room %','Callback',@obj.roomPercentChoose);
            set(obj.FigInfo.varChooseGrid,'ColumnSizes', [-1 -1],'RowSizes', [-1 -1 -1]);
            
            obj.FigInfo.varBChooseGrid = uiextras.Grid('Parent',obj.FigInfo.titleBar);
            obj.FigInfo.varTraj = uicontrol('Enable','on','Parent',obj.FigInfo.varBChooseGrid,'String','Traj','Callback',@obj.trajChoose);
            obj.FigInfo.varVisSpd = uicontrol('Enable','on','Parent',obj.FigInfo.varBChooseGrid,'String','Vis Spd','Callback',@obj.visSpdChoose);
            obj.FigInfo.varRunSpd = uicontrol('Enable','on','Parent',obj.FigInfo.varBChooseGrid,'String','Run Spd','Callback',@obj.runSpdChoose);
            obj.FigInfo.varDist = uicontrol('Enable','on','Parent',obj.FigInfo.varBChooseGrid,'String','Distance','Callback',@obj.distChoose);
            obj.FigInfo.varTotDist = uicontrol('Enable','on','Parent',obj.FigInfo.varBChooseGrid,'String','Total Distance','Callback',@obj.totDistChoose);
            obj.FigInfo.varRoomPercent = uicontrol('Enable','on','Parent',obj.FigInfo.varBChooseGrid,'String','Room %','Callback',@obj.roomPercentChoose);
            set(obj.FigInfo.varBChooseGrid,'ColumnSizes', [-1 -1],'RowSizes', [-1 -1 -1]);
            
            obj.FigInfo.lowContrast = axes('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.lowRoomLength = axes('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.lowGain = axes('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.noContrast = axes('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.normal = axes('Parent',obj.FigInfo.FigGrid);
            uiextras.Empty('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.highContrast = axes('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.highRoomLength = axes('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.highGain = axes('Parent',obj.FigInfo.FigGrid);
            
            set(obj.FigInfo.FigGrid,'ColumnSizes', [-1 -1 -1],'RowSizes', [-1 -1 -1]);
            makeAllPlots(obj)
        end
        function cellMinus(obj,~,~)
            if obj.ChosenCell>1
                obj.ChosenCell = obj.ChosenCell - 1;
                set(obj.FigInfo.cellDisplay,...
                'String',['Cell ' num2str(obj.ChosenCell) ': ' obj.Data.es.spikeIDs{obj.ChosenCell}]);
                makeAllPlots(obj)   
            end
        end
        function cellPlus(obj,~,~)
            if obj.ChosenCell<length(obj.Data.es.spikeIDs)
                obj.ChosenCell = obj.ChosenCell + 1;
                set(obj.FigInfo.cellDisplay,...
                'String',['Cell ' num2str(obj.ChosenCell) ': ' obj.Data.es.spikeIDs{obj.ChosenCell}]);
                makeAllPlots(obj)
            end
        end
        function makeAllPlots(obj)
            eval(['obj.variable = obj.Data.es.' obj.ChosenVar]);
            obj.plot_subplot(obj.Subset.normal, obj.FigInfo.normal);
            
            obj.plot_subplot(obj.Subset.noContrast, obj.FigInfo.noContrast);
            obj.plot_subplot(obj.Subset.lowContrast, obj.FigInfo.lowContrast);
            obj.plot_subplot(obj.Subset.highContrast, obj.FigInfo.highContrast);
            
            obj.plot_subplot(obj.Subset.lowGain, obj.FigInfo.lowGain);
            obj.plot_subplot(obj.Subset.highGain, obj.FigInfo.highGain);
            
            obj.plot_subplot(obj.Subset.lowRoomLength, obj.FigInfo.lowRoomLength);
            obj.plot_subplot(obj.Subset.highRoomLength, obj.FigInfo.highRoomLength);
            
            title(obj.FigInfo.normal,'NORMAL');
            title(obj.FigInfo.highRoomLength,'High Length');
            title(obj.FigInfo.noContrast,'Zero Contrast');
            title(obj.FigInfo.noContrast,'Low Contrast');
            title(obj.FigInfo.highContrast,'High Contrast');
            title(obj.FigInfo.highGain,'High Gain');
            title(obj.FigInfo.lowContrast,'Low Contrast');
            title(obj.FigInfo.lowGain,'Low Gain');
            title(obj.FigInfo.lowRoomLength,'Low Length');
            
        end
        function obj = plot_subplot(obj,t,plot_index)
            es = obj.Data.es;
            
            trialIDs = unique(es.trialID(t));
            for itr = 1:length(trialIDs)
                es.trialID(es.trialID==trialIDs(itr)) = itr;
            end
            
            plot(plot_index,es.trialID(t), obj.variable(t),'color',[.5 .5 .5])
            hold(plot_index,'on');
            plot(plot_index,es.trialID(es.lick & t), obj.variable(es.lick & t),'.','color',[.7 .7 .7])
            plot(plot_index,es.trialID(es.spikeTrain(:,obj.ChosenCell)>0 & t), obj.variable(es.spikeTrain(:,obj.ChosenCell)>0 & t),'ro')
            set(plot_index, 'box','off','TickDir','out','fontsize',14,'color','none')
            hold(plot_index,'off');
            axis tight
        end
        function obj = CalculateSubsets(obj)
            es = obj.Data.es;
            
            allCont = unique(es.contrast(~isnan(es.contrast)));
            if length(allCont)==4
                cont = allCont(3);
            elseif legnth(allCont)>1
                cont = allCont(end-1);
            else
                cont = allCont;
            end
            
            obj.Subset.all = es.traj>0 & es.smthBallSpd>5& ~isnan(es.smthBallSpd);
            obj.Subset.normal = obj.Subset.all & ...
                es.gain==1 & es.roomLength==1 & es.contrast==cont;
            obj.Subset.lowRoomLength = obj.Subset.all & ...
                es.gain==1 & es.roomLength<1 & es.contrast==cont;
            obj.Subset.highRoomLength = obj.Subset.all & ...
                es.gain==1 & es.roomLength>1 & es.contrast==cont;
            obj.Subset.lowGain = obj.Subset.all & ...
                es.gain<1 & es.roomLength==1 & es.contrast==cont;
            obj.Subset.highGain = obj.Subset.all & ...
                es.gain>1 & es.roomLength==1 & es.contrast==cont;
            obj.Subset.lowContrast = obj.Subset.all & ...
                es.gain==1 & es.roomLength==1 & es.contrast<cont & es.contrast>0;
            obj.Subset.highContrast = obj.Subset.all & ...
                es.gain==1 & es.roomLength==1 & es.contrast>cont;
            obj.Subset.noContrast = obj.Subset.all & ...
                es.contrast==0;
        end
        
        function trajChoose(obj,~,~)
            obj.ChosenVar = 'traj';
            set(obj.FigInfo.varTraj, 'Enable','off');
            set(obj.FigInfo.varVisSpd, 'Enable','on');
            set(obj.FigInfo.varRunSpd, 'Enable','on');
            set(obj.FigInfo.varDist, 'Enable','on');
            set(obj.FigInfo.varTotDist, 'Enable','on');
            set(obj.FigInfo.varRoomPercent, 'Enable','on');
            obj.makeAllPlots;
        end
        function visSpdChoose(obj,~,~)
            obj.ChosenVar = 'trajspeed';
            set(obj.FigInfo.varTraj, 'Enable','on');
            set(obj.FigInfo.varVisSpd, 'Enable','off');
            set(obj.FigInfo.varRunSpd, 'Enable','on');
            set(obj.FigInfo.varDist, 'Enable','on');
            set(obj.FigInfo.varTotDist, 'Enable','on');
            set(obj.FigInfo.varRoomPercent, 'Enable','on');
            obj.makeAllPlots;
        end
        function runSpdChoose(obj,~,~)
            obj.ChosenVar = 'ballspeed';
            set(obj.FigInfo.varTraj, 'Enable','on');
            set(obj.FigInfo.varVisSpd, 'Enable','on');
            set(obj.FigInfo.varRunSpd, 'Enable','off');
            set(obj.FigInfo.varDist, 'Enable','on');
            set(obj.FigInfo.varTotDist, 'Enable','on');
            set(obj.FigInfo.varRoomPercent, 'Enable','on');
            obj.makeAllPlots;
        end
        function distChoose(obj,~,~)
            obj.ChosenVar = 'distTrav';
            set(obj.FigInfo.varTraj, 'Enable','on');
            set(obj.FigInfo.varVisSpd, 'Enable','on');
            set(obj.FigInfo.varRunSpd, 'Enable','on');
            set(obj.FigInfo.varDist, 'Enable','off');
            set(obj.FigInfo.varTotDist, 'Enable','on');
            set(obj.FigInfo.varRoomPercent, 'Enable','on');
            obj.makeAllPlots;
        end
        function totDistChoose(obj,~,~)
            obj.ChosenVar = 'totDistTrav';
            set(obj.FigInfo.varTraj, 'Enable','on');
            set(obj.FigInfo.varVisSpd, 'Enable','on');
            set(obj.FigInfo.varRunSpd, 'Enable','on');
            set(obj.FigInfo.varDist, 'Enable','on');
            set(obj.FigInfo.varTotDist, 'Enable','off');
            set(obj.FigInfo.varRoomPercent, 'Enable','on');
            obj.makeAllPlots;
        end
        function roomPercentChoose(obj,~,~)
            obj.ChosenVar = 'trajPercent';
            set(obj.FigInfo.varTraj, 'Enable','on');
            set(obj.FigInfo.varVisSpd, 'Enable','on');
            set(obj.FigInfo.varRunSpd, 'Enable','on');
            set(obj.FigInfo.varDist, 'Enable','on');
            set(obj.FigInfo.varTotDist, 'Enable','on');
            set(obj.FigInfo.varRoomPercent, 'Enable','off');
            obj.makeAllPlots;
        end
    end
end
