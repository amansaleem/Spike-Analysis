function candidates = ProtocolFind (parnames, animallist, serieslist,options, xfilename)
% ProtocolFind finds experiments with specified active parameters
%
% L = ProtocolFind (parnames) searches the database specified in the global
% DIRS.data for experiments whose active parameters match all parameter names 
% provided in parnames. Place a tilde (~) in front of the parameter
% name to specify a 'passive' parameter, i.e. a parameter which should not
% be part of the active parameters list.
%
% L is a structure array with four fields: animal, iseries, iexp and description.
%
% L = ProtocolFind (parnames, animallist) lets you specify a list of
% animals. DEFAULT: [], which means all the animals in DIRS.data.
%
% L = ProtocolFind (parnames, animallist, serieslist) lets you specify a
% list of series. DEFAULT: [], which means all series.
%
% L = ProtocolFind (parnames, animallist, serieslist, options) to set options
% options(1): number of different values that active parameters must take (default: 2).
% options(2): when set to 1, excludes adaptation experiments (default:0).
%
%  L = ProtocolFind ({}, {}, [], [], xfilename) lets you specify the xfile
%  (do not add the .x at the end!). DEFAULT: '', which means all xfiles.
%
% Examples:
% L = ProtocolFind ({'c','diam'})
% L = ProtocolFind ({'~c','diam'})
% L = ProtocolFind({'seed'},{'Catz020'})
% L = ProtocolFind({'seed'},{'Catz020'},[3 5 7])
%
% part of the Spikes toolbox
%
% 2001 Matteo Carandini and Vincent Bonin
% 2001-11 MC added support for animallist, removed query on which animals to do
% 2002-01 VB added support for serieslist, ensures parnames/animallist are cellarrays
% 2002-02 VB added exclusion of certain active parameters
% 2002-05 MC added support for number of different values, made minor changes
% 2002-06 VB converted last ndiffvalues to an options vector. added adapt exclusion.
% 2003-02 VB corrected bug ListSeries was not called at every iteration
% 2003-10 MC safeguarded against ProtocolLoad fails, introduced report string
% 2009-11 MC improved text output
% 2011-03 AA excludes virtual reality data series from the list of series


reportstring = '';

defaultoptions = [2 0];

if nargin < 5
    xfilename = '';
end

if nargin < 4
    options = defaultoptions;
else
    nopt = length(options);
    options = [options defaultoptions(nopt+1:end)];
end

if nargin < 3
    serieslist = [];    
    if nargin < 2
        animallist = {};
    end
end

if ~iscell(parnames) || ~iscell(animallist)
    error('parameters ''parnames'' and ''animallist'' must be a cell array');
end

if ~isnumeric(serieslist)
    error('parameter ''serieslist'' must be a array');
end

global DIRS

if ~isfield(DIRS,'data')
    error('There must be a global DIRS, with a field ''data''');
end

% if ~isfield(DIRS,'spikes')
%     error('There must be a global DIRS, with a field ''spikes''');
% end

if ~isdir(DIRS.data)
    error(sprintf('Data directory %s does not exist', DIRS.data)); 
end

if isempty(animallist)
    strarrAnimals = ListAnimals(DIRS.data);
else
    strarrAnimals = animallist;  
end

if isempty(strarrAnimals)
    error(sprintf('Could not find any animal directories in %s\n',DIRS.data));
end

excludevec = strncmp(parnames,'~',1);
includelist = parnames(~excludevec);
excludelist = parnames( excludevec);

nanimals = length(strarrAnimals);

candidates = struct('animal','','iseries',[],'iexp',[],'description','');

% loop over animals
for ianimal = 1:nanimals
    animal = strarrAnimals{ianimal};
    disp(animal);
    if isempty(serieslist) % search all series
        strarrSeries = ListSeries(DIRS.data, animal);
    else % search specified series, if exist.
        strarrSeries = intersect(serieslist,ListSeries(DIRS.data, animal));
    end
    
    % AA: exclude VR series data (which are saved in a yymmdd format) from list of series
    strarrSeries(strarrSeries>100000)=[];
    
    
    % loop over series
    for iseries = strarrSeries(:)'
        disp(iseries);
        exptlist = ListExpts(DIRS.data, animal, iseries);
        
        % loop over experiments
        for iexp = exptlist(:)'
            failed = 0; 
            
            [ protocol, successflag ] = ProtocolLoad(animal,iseries,iexp); 
            if ~successflag || isempty(protocol)
                disp('WARNING: Could not read protocol');
                reportstring = sprintf('%s %s %d %d is not readable\n', reportstring, animal, iseries, iexp);
                failed = 1;
            end
            
            if ~failed && ~isempty(xfilename)
                if ~strcmp(protocol.xfile(1:end-2),xfilename)
                    failed = 1;
                end
            end
            
            if ~failed
                if protocol.adapt.flag && options(2) || isempty(protocol.activepars)
                    failed =1;
                end
            end
            
            if ~failed
                activepars = [protocol.activepars{:}]; % a flattened list
                % create array of active parameter names
                activeparnames = [];
                for iactivepar = 1:length(activepars)
                    activeparnames = [ activeparnames protocol.parnames(activepars(iactivepar)) ];
                end
                % loop over required parameters
                for iparname = 1:length(includelist)
                    parname = includelist(iparname);
                    iactivepar = strmatch(parname,activeparnames,'exact');
                    if isempty(iactivepar)
                        failed=1;
                    else
                        ipar = activepars(iactivepar);
                        if length(unique(protocol.pars(ipar,:))) < options(1)
                            failed = 1;
                        end
                    end
                end
                % loop over exclusion list
                for iparname = 1:length(excludelist)
                    parname = excludelist{iparname};
                    % removing tilde
                    parname = parname(2:end);
                    if ismember(parname,activeparnames)
                        failed=1;
                    end
                end
            end
            
            if ~failed
                % this is a hit!
                candidates(end+1).animal = animal;
                candidates(end).iseries = iseries;
                candidates(end).iexp = iexp;
                candidates(end).description = protocol.description;
                disp('----------------------> Match found!');
                reportstring = sprintf('%s %s %d %d -- %s\n', ...
                    reportstring, animal, iseries, iexp, protocol.description);
                
            end

        end % for iexp
    end % for iseries
end % for ianimal

% first record is empty
candidates=candidates(2:end);

fprintf(1,reportstring);