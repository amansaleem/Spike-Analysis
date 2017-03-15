classdef ExperimentTag
    % ExperimentTag specifies an animal, experiment, and series
    %
    % Properties:
    %   animal  -- must be a string
    %   iseries -- must be an integer positive number
    %   iexp    -- must be an integer positive number
    %
    % Methods:
    %   ShortName gives a string that combines animal, iseries, iexp
    %   LongName gives a longer string with more information
    % 
    % 2012-01 MC
    
    properties
        animal
        iseries
        iexp
    end
    
    methods
        
        function ET = ExperimentTag(animal,iseries,iexp)     
            
            if ~ischar(animal)
                error('animal name must be a string');
            end
            if ~isnumeric(iseries)||(iseries~=round(iseries))||(iseries<=0)
                error('series must be an integer positive number');
            end
             if ~isnumeric(iexp)||(iexp~=round(iexp))||(iexp<=0)
                error('experiment must be an integer positive number');
            end
           
            ET.animal = animal;
            ET.iseries = iseries;
            ET.iexp = iexp;            
        end
        
        function StringTag = ShortName(ET)
            StringTag = sprintf('%s-%d-%d',ET.animal,ET.iseries,ET.iexp);
        end
        
        function StringTag = LongName(ET)
            StringTag = sprintf('Animal %s series %d experiment %d',ET.animal,ET.iseries,ET.iexp);
        end
        
    end
    
end

