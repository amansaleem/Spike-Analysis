function strarrAnimals = ListAnimals(datadir)
% ListAnimals lists all the animals in a data directory
% 
% strarrAnimals = ListAnimals gives the list of animals in DIRS.data. The
% list is a cell array of strings.
%
% strarrAnimals = ListAnimals(datadir) lets you specify the directory.
%
% part of Spikes toolbox

% 2001-03 Matteo Carandini
% 2002-04 MC sorted the output
% 2010-02 MC added support for global DIRS

global DIRS

if nargin<1
    datadir = [];
end

if isempty(datadir) && isfield(DIRS,'data')
    datadir = DIRS.data;
end

if ~exist(datadir,'dir')
   error(sprintf('Data directory does not exist: \n%s',datadir));
end

dd = dir(datadir);

strarrAnimals = {dd([dd.isdir] == 1).name};
strarrAnimals = setdiff( strarrAnimals, {'.', '..','pfiles','xfiles'});

ndirs = length(strarrAnimals);
checked = zeros(ndirs,1);
for idir = 1:ndirs
   % check that in the dir there is a log file
   dirname = strarrAnimals{idir};
   if exist(fullfile(datadir,dirname,[ dirname '.txt']),'file')
      checked(idir) = 1;
   end
end

strarrAnimals = {strarrAnimals{find(checked)}}; %#ok<FNDSB>

[strarrUpperAnimals, ii] = sort(upper(strarrAnimals));

strarrAnimals = strarrAnimals(ii);
