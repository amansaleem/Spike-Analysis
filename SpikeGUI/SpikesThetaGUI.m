classdef SpikesThetaGUI < SpikeGUIController
    
    properties
        
    end
    
    methods
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
                obj.data.es = VRLoadMultipleExpts(obj.animal, obj.iseries, obj.exptList,'SPIKES_THETA');
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
        
    end 
end