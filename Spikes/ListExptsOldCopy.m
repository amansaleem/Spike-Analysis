function Exptlist = ListExpts(datadir, stranimal, iseries)
% ListExpts lists all the expts that were run for a given series
% 
% Exptlist = ListExpts(stranimal, iseries)
% returns a vector
% 
% Exptlist = ListExpts(datadir, stranimal, iseries)
% returns a vector
%
% 2001-02 Matteo Carandini
% 2007-08 MC added case of two inputs
% 
% part of Spikes

global DIRS

if nargin == 2

    iseries = stranimal;
    stranimal = datadir;
    datadir = DIRS.data;
end

if ~exist(datadir,'dir')
   error('Data directory does not exist');
end

if ~isstr(stranimal)
   error('Animal must be a string');
end

strseries = num2str(iseries);

Exptdir = dir(fullfile(datadir, stranimal, strseries));

if isempty(Exptdir)
   return
end

% if we are looking for the actual data:
strarrExpt = {Exptdir(find([Exptdir.isdir] == 1)).name};
strarrExpt = setdiff( strarrExpt, {'.', '..'});
if isempty(strarrExpt), 
   % set(pickerfig,'pointer','arrow');
   return; 
end
ne = length(strarrExpt);
% convert strarrExpt to a matrix to sort and then convert back!
for ie = 1:ne
   Exptlist(ie) = sscanf(strarrExpt{ie},'%d');
end
Exptlist = sort(Exptlist);
