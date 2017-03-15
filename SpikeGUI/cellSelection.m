classdef cellSelection < handle & hgsetget
    
    properties
        ChosenCell = 1;
        NumCells   = 1;
        CellList   = 1;
        CellListString = 'load';
        CellListFull=1;
        NumCellsFull=1;
    end
    
    methods
        function obj = cellSelection
        end
        function obj = cellMinus(obj)
            if obj.ChosenCell>1
                obj.ChosenCell = obj.ChosenCell - 1;
            end
        end
        function obj = cellPlus(obj)
            if obj.ChosenCell<obj.NumCells
                obj.ChosenCell = obj.ChosenCell + 1;
            end
        end
        function obj = defineCellInfo(obj, spikeIDs)
            obj.NumCells = length(spikeIDs);
            obj.CellList = 1:obj.NumCells;
            obj.CellListFull = 1:obj.NumCells;
            obj.NumCellsFull = length(spikeIDs);
            obj.CellListString = [];
            for icell = 1:obj.NumCells
                obj.CellListString{icell} = [num2str(icell) ': ' spikeIDs{icell}] ;
            end
        end
        function obj = redefineCellInfo(obj, newCellList)
            obj.NumCells = length(newCellList);
            obj.CellList = newCellList;
            oldListString = obj.CellListString;
            obj.CellListString = [];
            for icell = 1:obj.NumCells
                obj.CellListString{icell} = [num2str(newCellList(icell)) ': ' oldListString{newCellList(icell)}];
            end
        end
    end
end