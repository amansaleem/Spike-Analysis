classdef SpikeDisplayGUI < handle
    
    properties
        FigInfo
        Data
        CellInfo 
        ChosenCell = 1;
        ChosenVar  = 'traj';
        ChosenVarX  = 'trial';
        Subset
    end
    
    methods
        function obj = SpikeDisplayGUI
        end
    end
    methods (Abstract)
        makeAllPlots(obj,icell)
        createUI
    end
end