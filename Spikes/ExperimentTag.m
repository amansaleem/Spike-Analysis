classdef ExperimentTag
    % ExperimentTag specifies an animal, experiment, and series
    %
    % myExptTag = ExperimentTag(animal,iseries,iexp)
    % myExptTag = ExperimentTag('aaaa-ss-ii')
    % myExptTag = ExperimentTag(struct), where struct has fields animal, iseries, iexp
    %
    % Properties:
    %   animal  -- must be a string
    %   iseries -- must be an integer positive number
    %   iexp    -- must be an integer positive number
    %
    % Methods:
    %   ShortName gives a string 'aaaa-ss-ii' combining animal, iseries, iexp
    %   LongName gives a longer string with more information
    %   
    %
    % 2012-01 MC
    % 2012-02 MC added option to read string
    % 2012-02 MC added option to read struct
    
    properties
        animal
        iseries
        iexp
    end
    
    
     methods(Static)
        
        function ET = ExperimentTag(varargin)
            
            switch nargin
                case 3
                    TheAnimal  = varargin{1};
                    TheSeries  = varargin{2};
                    TheExp     = varargin{3};
                case 1
                    v = varargin{1};
                    if ischar(v)
                        TheFields  = textscan(v,'%s %d %d','delimiter','-');
                        TheAnimal  = TheFields{1}{1};
                        TheSeries  = TheFields{2};
                        TheExp     = TheFields{3};
                    elseif isstruct(v) && isfield(v,'animal') && isfield(v,'iseries') && isfield(v,'iexp')
                        TheAnimal  = v.animal;
                        TheSeries  = v.iseries;
                        TheExp     = v.iexp;
                    else
                        error('Do not understand inputs');
                    end
                    
                otherwise
                    error('Do not understand inputs');
                    
            end
            
            if ~ischar(TheAnimal)
                disp(TheAnimal)
                error('animal name must be a string');
            end
            if ~isnumeric(TheSeries)||(TheSeries~=round(TheSeries))||(TheSeries<=0)
                error('series must be an integer positive number');
            end
            if ~isnumeric(TheExp)||(TheExp~=round(TheExp))||(TheExp<=0)
                error('experiment must be an integer positive number');
            end
            
            ET.animal = TheAnimal;
            ET.iseries = TheSeries;
            ET.iexp = TheExp;
        end
        
       
     end
     
     methods
        
    
        function StringTag = ShortName(ET)
            StringTag = sprintf('%s-%d-%d',ET.animal,ET.iseries,ET.iexp);
        end
        
        function StringTag = LongName(ET)
            StringTag = sprintf('Animal %s series %d experiment %d',ET.animal,ET.iseries,ET.iexp);
        end
        
    end
    
end

