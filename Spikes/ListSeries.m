function serieslist = ListSeries(datadir, animal)
% ListSeries lists all the series that were run for a given animal
% 
% serieslist = ListSeries(animal)
% returns a vector
%
% serieslist = ListSeries(datadir, animal)
% returns a vector
%
% 2001-02 Matteo Carandini
% 2007-08 MC added case of one input
% 
% part of Spikes

global DIRS
    
if isempty(DIRS)
    error('Must define DIRS -- perhaps you should run SetDefaultDirs');
end

if nargin == 1

    animal = datadir;
    datadir = DIRS.data;
end

if ~exist(datadir,'dir')
   error('Data directory does not exist');
end

if ~ischar(animal)
   error('Animal must be a string');
end

seriesdir = dir(fullfile(datadir,animal));

strarrSeries = {seriesdir(find([seriesdir.isdir] == 1)).name};
strarrSeries = setdiff(strarrSeries,{'.','..'});

ns = length(strarrSeries);

% convert strarrSeries to a matrix to sort and then convert back!
serieslist = [];

if ns > 0
    for is = 1:ns
        thisentry =  sscanf(strarrSeries{is},'%d');
        if ~isempty(thisentry), serieslist(is) = thisentry; end
    end
    if ~isempty(serieslist)
        serieslist = sort(serieslist);
    end
end
