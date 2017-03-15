function hax = UnitPlotCycles( unit, resolution, stims, protocol)
% UnitPlotCycles plots cycle hystograms
%
% UnitPlotCycles( unit ) plots the cycle histograms for all stimuli, arranged 
% in two columns. 
%
% UnitPlotCycles( unit, resolution ) lets you specify the resolution in seconds
% (default: 0.01). Can be also a vector containing different resolutions for 
% different stimuli
% 
% hax = UnitPlotCycles( unit, resolution, stims)lets you specify a matrix of stimuli
% (in case you don't want all stimuli or you want to arrange them in a particular way).
% Entries with 0 indicate an empty histogram.
%
% hax = UnitPlotCycles( unit, resolution, stims, protocol) lets you specify a 
% protocol (can be a model protocol, doesn't have to be the protocol of a real
% experiment
%
% hax = UnitPlotCycles(  ) returns a matrix of pointers to the axes.
%
% part of Spikes
% 2001-02 MC
% 2001-08 VM added screenspec, which is an input to unitgetcycles
% 2001-10-10 VM	removed screenspec, now uses protocol.estfreqs
%						ProtocolLoad called with (...,'loadscreen')
% 2001-04-20 VM	small changes, added protocol as an input, let resolution be a vector 
% 2003-03 VM made it fit for traces

% 
% -------------- Example:
%
% global DIRS
% DIRS.spikes = 'f:/Spikes';
% DIRS.data = 'f:/Cat';
%
% % --- Load a spatial frequency experiment (period is the same for all):
% unit = UnitLoad( DIRS.spikes, 'CATZ008', 5, 1, 1, 2);
% figure; 
% UnitPlotCycles(unit, 0.005);
% 
% % --- Load a temporal frequency experiment (period varies from stim to stim):
% unit = UnitLoad( DIRS.spikes, 'CATZ008', 5, 2, 1, 2);
% figure; 
% UnitPlotCycles(unit, 0.005);
%
% % --- Plot a fancy version with specific stimuli in specific positions:
% figure; 
% UnitPlotCycles(unit, 0.005, [7 0 0; 3 2 1; 3 3 3]);

% Decide what datatype to use
datatype = unit.datatype;

% See if datatype is right and data exist
if ~strcmp(unit.datatype,'spiketimes') || isempty(unit.spiketimes)
    warning('UnitPlotCycles:NotImplemented','Can''t plot cycles because there are no spikes');
    return
end



animal 	= unit.animal;
iseries 	= unit.iseries;
iexp 		= unit.iexp;
ichan 	= unit.ichan;
icell		= unit.icell;

if nargin < 4 || isempty(protocol)
   protocol = ProtocolLoad( animal, iseries, iexp, 'loadscreen');
end

data = getfield(unit,datatype);
nstim = size(data,1);
nrpts	= size(data,2);

if nargin < 3 || isempty(stims)
   nrows = ceil(nstim/3);
   ncols = 3;
   stims = zeros(nrows, ncols);
   stims(1:nstim) = 1:nstim;
   labelflag = 'labels';
else
   if ndims(stims)~=2
      error('Stims must be a  2-D matrix');
   end
   if any(stims>nstim | stims<0)
      error('Stimulus matrix contains invalid stimulus numbers');
   end
   [nrows,ncols] = size(stims);
   labelflag = 'nolabels';
end

if nargin < 2 || isempty(resolution)
   resolution = 0.01; % s
else
   if any(resolution < 0.001) || any(resolution > 0.2)
      error('Resolution cannot be less than 0.002 s or more than 0.2 seconds (for arbitrary reasons)');
   end
end

%% Get the cycle data

[rrcell periods] = UnitGetCycles( unit, protocol, resolution, [0]); 

%% Graphics

figure;
hax = zeros(nrows,ncols); % histo axes

for irow = 1:nrows
   for icol = 1:ncols
      istim = stims(irow,icol);
      
      hax(irow,icol) = subplot(nrows,ncols,ncols*(irow-1)+icol);
      if istim<=0
         set(hax(irow,icol),'visible','off');
      else
         if strcmp(labelflag,'labels')
            text(0,0.5,num2str(istim),...
               'hori','right','vert','middle','units','norm','fontsize',18);
         end
         hold on;
         nbins = length(rrcell{istim});
         fillplot(linspace(0,1,nbins), zeros(1,nbins), rrcell{istim},[.5 .5 .5]);
         set(hax(irow,icol),'xtick',1,'xticklabel',num2str(periods(istim),2));
      end
   end
end

set(hax,'xlim',[-inf inf],'ylim',[0 inf],'ytick',[]);
set(hax,'box','off','xcolor','k','ycolor','w')
[mn, mx] = matchy(hax);

lasthax = find(stims(:)>0, 1, 'last' );

set(hax(lasthax),...
   'ycolor','k','ytick',[floor(mx)],'yaxislocation','right');

if length(unique(periods))==1
   set(hax,'xtick',[]);
   set(hax(lasthax),'xtick',1,'xticklabel',num2str(unique(periods),2));
end

axes(hax(1,round(ncols/2)));
title(['Unit ' animal '.' num2str(ichan) '.' num2str(icell) '   Expt ' num2str(iseries) '.' num2str(iexp) ]);

[mn, mx] = matchy(hax);
axes(hax(1));
text(1,1,num2str(round(mx)),'units','norm','hori','left','vert','middle');
text(1,0,num2str((periods(1))),'units','norm','hori','left','vert','middle');

set(hax,'plotboxaspectratio',[1 1 1]);

%% the description

desc = sprintf('%s Chan %d Cell %d Series %d Exp %d ',unit.animal,unit.ichan,unit.icell, unit.iseries, unit.iexp);
set(gcf,'menubar','none','numbertitle','off','name',desc);
if nargin < 3 || strcmp(labelflag,'labels')
   set(gcf,'position',[100   100   250   500]);
end

