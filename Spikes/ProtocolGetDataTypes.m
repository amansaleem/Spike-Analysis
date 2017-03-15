function DataTypes = ProtocolGetDataTypes(p)
% checks out what kind of data were acquired during an experiment
%
% DataTypes = ProtocolGetDataTypes(protocol)
% returns
% DataTypes.Traces          = true or false;
% DataTypes.TDT             = true or false;
% DataTypes.CerebusSnippets = true or false;
% DataTypes.CerebusTraces   = true or false;
%
% 2010-03 Matteo
% 2010-06 MC if there is a Nev but no Ns3, then ignore it (photodiode only)
% 2010-10 ND had the check for blackrock traces look to see if there are
% appropriately named files in the series and experiment subfolders when saving by repeat
% 2010-12 MS made 'BlackRock snippets' button light up whenever ns3, ns4, or ns5 file is recorded


global DIRS

%% Start pessimistic

DataTypes.Traces          = false;
DataTypes.TDT             = false;
DataTypes.CerebusSnippets = false;
DataTypes.CerebusTraces   = false;

%% see if there are Multispike or Traces data

thedir = fullfile(DIRS.data,p.animal,int2str(p.iseries),int2str(p.iexp));
thefile = sprintf('%s_%d_%d_*.mat',p.animal, p.iseries, p.iexp);
dd = dir(fullfile(thedir,thefile));
DataFileExists = 0;
for id = 1:length(dd)
    if isempty(strfind(dd(id).name,'Michigan'))
        DataFileExists = 1;
        break;
    end
end
if DataFileExists
    foo = load(fullfile(thedir,dd(id).name),'stimresps');
    if isfield(foo.stimresps,'data')
       DataTypes.Traces = true;
    end
end

%% see if there are TDT data

thedir = fullfile(DIRS.data,p.animal,int2str(p.iseries),int2str(p.iexp));
thefile = sprintf('%s_%d_%d_*-Michigan.mat',p.animal, p.iseries, p.iexp);
dd = dir(fullfile(thedir,thefile));
if numel(dd)>0
    DataTypes.TDT = true;
end

%% see if there are Blackrock snippet data

thedir = fullfile(DIRS.Cerebus,p.animal);
TheNevFile = sprintf('u%03.f_%03.f.nev',p.iseries,p.iexp);
if exist(sprintf('u%03.f_%03.f.ns3',p.iseries,p.iexp),'file')
    TheNsFile = sprintf('u%03.f_%03.f.ns3',p.iseries,p.iexp);
elseif exist(sprintf('u%03.f_%03.f.ns4',p.iseries,p.iexp),'file')
    TheNsFile = sprintf('u%03.f_%03.f.ns4',p.iseries,p.iexp);
elseif exist(sprintf('u%03.f_%03.f.ns5',p.iseries,p.iexp),'file')
    TheNsFile = sprintf('u%03.f_%03.f.ns5',p.iseries,p.iexp);
else
    TheNsFile = sprintf('xxx');
end

% if there is a Nev but no Ns*, then the Nev just has the photodiode times
ddNev = dir(fullfile(thedir,TheNevFile));
ddNs  = dir(fullfile(thedir,TheNsFile));
if numel(ddNev)>0 && numel(ddNs)>0
    DataTypes.CerebusSnippets = true;
end

%% see if there are Blackrock traces data

thedir = fullfile(DIRS.Cerebus,p.animal);
thefile = sprintf('u%03.f_%03.f.ns5',p.iseries,p.iexp);
thesmallerfile = sprintf('u%03.f_%03.f.ns4',p.iseries,p.iexp); % added by ND 20101019
% put in here check for the subfolders - ND201021
thesubdir = fullfile(DIRS.Cerebus,p.animal,num2str(p.iseries),num2str(p.iexp));
therptfiles = [BuildFileName(thesubdir,p.animal,p.iseries,p.iexp) '*.ns5'];
thesmallerrptfiles = [BuildFileName(thesubdir,p.animal,p.iseries,p.iexp) '*.ns4'];

CerebusTracesByExpt = false;
dd = dir(fullfile(thedir,thefile));
if numel(dd)>0
    CerebusTracesByExpt = true;
else % added by ND 20101019
    ee = dir(fullfile(thedir,thesmallerfile));
    if numel(ee)>0
        CerebusTracesByExpt = true;
    end
end
% added by ND 201021
CerebusTracesByRpt = false;
ff = dir(therptfiles);
if numel(ff)>0
    CerebusTracesByRpt = true;
else
    gg = dir(thesmallerrptfiles);
    if numel(gg)>0
        CerebusTracesByRpt = true;
    end
end

DataTypes.CerebusTraces = CerebusTracesByExpt | CerebusTracesByRpt;

