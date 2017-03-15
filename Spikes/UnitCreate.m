function u = UnitCreate( pseudounit )
% UnitCreate creates an empty unit structure
%
% u = UnitCreate
%
% u = UnitCreate( pseudounit ) creates the unit structure using fields from
% the possibly malformed pseudounit
%
% 2008-03 Matteo Carandini
% 2008-03 LB checked number of arguments
% 2010-06 added field traces
% 2011-03 MS added fields prestimdurs, prespiketimes, author, sortMethod, sortDate (described in google document 'The Unit Structure')
% 2011-03 MS added fields waveformAll, waveformMean, waveformStd, isolDist, LRatios (described in google document 'The Unit Structure')
% 2013-08 MC and ND added fields poststimdurs and postspiketimes

if nargin < 1
    pseudounit = [];
end

% AZ20090710
% u = struct(                     ...
%            'animal'      ,{'' },...
%            'ichan'       ,{NaN},...
%            'icell'       ,{NaN},...
%            'iseries'     ,{NaN},...
%            'iexp'        ,{NaN},...
%            'stimdurs'    ,{   },...
%            'timestamp'   ,{NaN},...
%            'prototype'   ,{   },...
%            'neighborhood',{   },...
%            'spiketimes'  ,{   },...
%            'datatype'    ,{'' },...
%            'source'      ,{'' } ...
%                                    );

u.animal        = '';
u.iseries       = NaN; 
u.iexp          = NaN; 
u.ichan         = NaN;
u.icell         = NaN; 
u.id            = '';
u.datatype      = ''; % 'traces' or 'spiketimes'
u.traces        = []; 
u.timestamp     = NaN; 
u.nstim         = [];
u.nrepeats      = [];
u.stimdurs      = []; 
u.prestimdurs   = []; 
u.spiketimes    = [];
u.prespiketimes = []; 
u.prototype     = []; 
u.sampledur     = [];
u.waveformAll   = [];
u.waveformMean  = [];
u.waveformStd   = [];
u.isolDist      = [];
u.lRatios       = [];
u.neighborhood  = [];
u.source        = '';
u.author        = ''; 
u.sortMethod    = ''; 
u.sortDate      = ''; 

% MC and ND added these two fields on 2013-08-29
u.poststimdurs   = [];
u.postspiketimes = [];

if ~isempty(pseudounit)
    names = fieldnames(u);
    for iname = 1 : length(names)
        if isfield(pseudounit, names{iname})
            u.(names{iname}) = pseudounit.(names{iname});
        end
    end
end