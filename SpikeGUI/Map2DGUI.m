classdef Map2DGUI < SpikeDisplayGUI & hgsetget
    
    properties
        variable
        variableB
        maps2d
%         showErrorBars = 1;
    end
    
    methods
        % add an update button, this can be called by central
        function obj = Map2DGUI(data)
            obj.Data.es = data;
            obj = obj.calculateSubsets;
            obj = chooseVarUI(obj);
            
            obj = obj.calculateMaps;
            
            obj = obj.createUI;
        end
        function obj = createUI(obj)
            obj.FigInfo.fig = figure('Name', 'Map 2D Display',...
                'Position', [700 100 500 900]);
            obj.FigInfo.layout = uiextras.VBox('Parent', obj.FigInfo.fig);
            obj.FigInfo.titleBar = uiextras.HBox('Parent', obj.FigInfo.layout);
            obj.FigInfo.FigGrid = uiextras.Grid('Parent', obj.FigInfo.layout);
            obj.FigInfo.layout.Sizes = [80 -1];
            set(obj.FigInfo.titleBar, 'Spacing',20)
            
            obj.FigInfo.cellSection = uiextras.VBox('Parent', obj.FigInfo.titleBar);
            obj.FigInfo.cellDisplay = uicontrol('Style','text','Parent', obj.FigInfo.cellSection,'fontsize',12,...
                'String',['Cell ' num2str(obj.ChosenCell) ': ' obj.Data.es.spikeIDs{obj.ChosenCell}]);
            uiextras.Empty('Parent', obj.FigInfo.cellSection);
            
            obj.FigInfo.normal = axes('Parent',obj.FigInfo.FigGrid);
            uiextras.Empty('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.noContrast = axes('Parent',obj.FigInfo.FigGrid);
            
            obj.FigInfo.lowContrast = axes('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.lowRoomLength = axes('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.lowGain = axes('Parent',obj.FigInfo.FigGrid);
            
            obj.FigInfo.highContrast = axes('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.highRoomLength = axes('Parent',obj.FigInfo.FigGrid);
            obj.FigInfo.highGain = axes('Parent',obj.FigInfo.FigGrid);
            
            set(obj.FigInfo.FigGrid,'ColumnSizes', [-1 -1 -1],'RowSizes', [-1 -1 -1]);
            
        end
        % Move figure related components to createUI
        function obj = chooseVarUI(obj)
            
            VarTypes = {'Traj','Distance','Total Distance','Room %','Vis Spd','Run Spd'};
            [s,v] = listdlg('PromptString','Select the X variable:',...
                'SelectionMode','single',...
                'ListString',VarTypes);
            if v
                switch s
                    case 1; obj.ChosenVarX = 'traj';
                    case 2; obj.ChosenVarX = 'distTrav';
                    case 3; obj.ChosenVarX = 'totDistTrav';
                    case 4; obj.ChosenVarX = 'trajPercent';
                    case 5; obj.ChosenVarX = 'trajspeed';
                    case 6; obj.ChosenVarX = 'smthBallSpd';
                end
            else
                obj.ChosenVarX = 'trajPercent';
            end
            
            VarTypes = {'Traj','Distance','Total Distance','Room %','Vis Spd','Run Spd'};
            [s,v] = listdlg('PromptString','Select the Y variable:',...
                'SelectionMode','single',...
                'ListString',VarTypes);
            if v
                switch s
                    case 1; obj.ChosenVar = 'traj';
                    case 2; obj.ChosenVar = 'distTrav';
                    case 3; obj.ChosenVar = 'totDistTrav';
                    case 4; obj.ChosenVar = 'trajPercent';
                    case 5; obj.ChosenVar = 'trajspeed';
                    case 6; obj.ChosenVar = 'smthBallSpd';
                end
            else
                obj.ChosenVar = 'smthBallSpd';
            end
            
        end
        
        function obj = calculateMaps(obj)
            obj.maps2d = [];
            eval(['obj.variable(:,1) = obj.Data.es.' obj.ChosenVarX ';']);
            eval(['obj.variable(:,2) = obj.Data.es.' obj.ChosenVar ';']);
            obj = obj.getTwoMap('all');
            obj = obj.getTwoMap('normal');
            obj = obj.getTwoMap('lowRoomLength');
            obj = obj.getTwoMap('highRoomLength');
            obj = obj.getTwoMap('lowGain');
            obj = obj.getTwoMap('highGain');
            obj = obj.getTwoMap('noContrast');
            obj = obj.getTwoMap('lowContrast');
            obj = obj.getTwoMap('highContrast');
            
        end
        function obj = getTwoMap(obj,type)
            progressMessage = waitbar(0,['Processing 2D maps for ' type '...']);
            eval(['obj.maps2d.' type ' = twoDimMap;']);
            eval(['datapts = sum(obj.Subset.' type ');']);
            if datapts>100
                %This is a random number for now
                eval(['obj.maps2d.' type ' = obj.maps2d.' type '.trainSpikeMap( obj.variable(obj.Subset.' type ',:), obj.Data.es.spikeTrain(obj.Subset.' type ',:) );']);
            end
            close(progressMessage)
        end
        function makeAllPlots(obj,icell)
%             icell = obj.ChosenCell;
%             display(['Chosen cell: ' num2str(icell)]);
            eval(['obj.variable = obj.Data.es.' obj.ChosenVar ';']);
            
            set(obj.FigInfo.cellDisplay,...
                'String',['Cell ' num2str(icell) ': ' obj.Data.es.spikeIDs{icell}]);
            
            plot2Dmap(obj.maps2d.normal,icell,obj.FigInfo.normal)
%             axis off
            title(obj.FigInfo.normal,'NORMAL');
            
            plot2Dmap(obj.maps2d.noContrast,icell,obj.FigInfo.noContrast)
%             axis off
            title(obj.FigInfo.noContrast, 'Zero Contrast');
            
            plot2Dmap(obj.maps2d.lowRoomLength,icell,obj.FigInfo.lowRoomLength)
            title(obj.FigInfo.lowRoomLength, 'Low RoomLength');
            plot2Dmap(obj.maps2d.highRoomLength,icell,obj.FigInfo.highRoomLength)
            title(obj.FigInfo.highRoomLength, 'High RoomLength');
            
            plot2Dmap(obj.maps2d.lowContrast,icell,  obj.FigInfo.lowContrast )
            title(obj.FigInfo.lowContrast, 'Low Contrast');
            plot2Dmap(obj.maps2d.highContrast,icell,  obj.FigInfo.highContrast )
            title(obj.FigInfo.highContrast, 'High Contrast');
            
            plot2Dmap(obj.maps2d.lowGain,icell,obj.FigInfo.lowGain)
            title(obj.FigInfo.lowGain, 'High Gain');
            plot2Dmap(obj.maps2d.highGain,icell,obj.FigInfo.highGain)
            title(obj.FigInfo.highGain, 'High Gain');
            
            
            
        end
        function obj = plot_subplot(obj,plot_index,icell)
            plot2Dmap(obj.model,icell,plot_index)
        end
        function obj = calculateSubsets(obj)
            es = obj.Data.es;
            
            allCont = unique(es.contrast(~isnan(es.contrast)));
            if length(allCont)==4
                cont = allCont(3);
            elseif length(allCont)>1
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
    end
end
