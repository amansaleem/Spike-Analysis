function Exptlist = ListExpts(datadir, stranimal, iseries, datatype)
% ListExpts lists all the expts that were run for a given series
%
% Exptlist = ListExpts(datadir, stranimal, iseries, 'TDT')
% lists all TDT/michigan experimnents
%
% Exptlist = ListExpts(datadir, stranimal, iseries, 'BR')
% lists all BlackRock experimnents
%
% Exptlist = ListExpts(stranimal, iseries)
% returns a vector
%
% Exptlist = ListExpts(datadir, stranimal, iseries)
% returns a vector
%
% 2001-02 Matteo Carandini
% 2007-08 MC added case of two inputs
% 2010-03 AS added datatype to list Expts from either TDT or BlackRock
%
% part of Spikes

global DIRS

if nargin == 2
    
    iseries = stranimal;
    stranimal = datadir;
    datadir = DIRS.data;
end

% AS: if we want to list expts of certain type TDT / BlackRock
if (nargin == 4 | (nargin == 3 & isstr(iseries)))
    if nargin == 3
        datatype = iseries;
        iseries = stranimal;
        stranimal = datadir;
        datadir = DIRS.data;
    end
    switch datatype
        case 'TDT'
            tdt_flag = 1;
            br_flag  = 0; %BlackRock
        case 'BR'
            tdt_flag = 0;
            br_flag  = 1; %BlackRock
    end
else
    tdt_flag = 0;
    br_flag  = 0; %BlackRock
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


% AS: Find TDT expts
if tdt_flag
    new_Exptlist = [];
    for exptID = 1:length(Exptlist)
        thedir = fullfile(datadir, stranimal, strseries, int2str(Exptlist(exptID)));
        thefile = sprintf('%s_%s_%d_*-Michigan.mat',stranimal, strseries, Exptlist(exptID));
        dd = dir(fullfile(thedir, thefile));
        if numel(dd)>0
            new_Exptlist = [new_Exptlist Exptlist(exptID)];
        end
    end
    Exptlist = new_Exptlist;
end

% AS: Find BlackRock expts
if br_flag
    new_Exptlist = [];
    for exptID = 1:length(Exptlist)
        thedir = fullfile(DIRS.Cerebus, stranimal);
        thefile = sprintf('u%03.f_%03.f.nev',str2num(strseries), Exptlist(exptID));
        dd = dir(fullfile(thedir,thefile));
        if numel(dd)>0
            new_Exptlist = [new_Exptlist Exptlist(exptID)];
        end
    end
    Exptlist = new_Exptlist;
end