function hax = UnitPlotHistos( unit, stimlist, stdflag )
% UnitPlotHistos plots the histograms of responses in a Unit data structure
%
% UnitPlotHistos( unit )
%
% UnitPlotHistos( unit, stimlist )
%
% UnitPlotHistos( unit, stimlist, stdflag )
%
% part of Spikes
%
% 2000-11 MC
% 2001-10 VM plots also errorbars
% 2001-12 MC changed to active filtering, changed plotting
% 2002-12 MC changed to output the axes
% 2003-03 VM small changes to make it plot also traces
% 2003-03 VM made it fit for traces
% 2003-07 MC allowed for stimlist in inputs, added stdflag

% unit = UnitLoad( DIRS.spikes, 'CATZ008', 5, 9, 1, 2);

if nargin == 3
    stdflag = 1;
else
    stdflag = 0;
end

% See if datatype is right and data exist
if ~strcmp(unit.datatype,'spiketimes') || isempty(unit.spiketimes)
    warning('UnitPlotHistos:NotImplemented','Can''t plot histograms because there are no spikes');
    return
end

[R, E, resolution] = UnitGetRates( unit, 0.015, 'active' ); % was 0.03 MC for catz028

fprintf(1,'Obtained firing rates with resolution %2.1f ms\n',1000*resolution);



nstim = length(R);

if nargin<2
    nax = nstim;
    stimlist = 1:nstim;
    ncols = ceil(nstim/20);
else
    nax = length(stimlist);
    ncols = 1;
end

figure; clf; hax = zeros(nax,1); % histo axes

nrows = ceil(nax/ncols);

for iax = 1:nax
    
    istim = stimlist(iax);
    
    hax(iax) = subplot(nrows,ncols,iax);
    text(0,0.5,num2str(istim),'hori','right','vert','middle','units','norm','fontsize',18);
    hold on;
    
    n = length(R{istim});
    yy = [R{istim};R{istim}];
    yy = yy(:);
    ee = [E{istim};E{istim}];
    ee = ee(:);
    if stdflag
        % turn the error info from sem into std
        ee = ee*sqrt(unit.nrepeats);
    end
    xx = [1:n; 1:n];
    xx = xx(:);
    xx = [ 0; xx(1:end-1) ];
    
    plot(xx,yy+ee,'color',[.6 .6 .6]); hold on
    if stdflag
        fillplot(xx,zeros(size(xx)),yy+ee,'r');
    end
    fillplot(xx,zeros(size(xx)),yy,[0 0 0]);
    
    maxdur = max(unit.stimdurs(istim,:));
    text(1,0,num2str(round(maxdur*100)/100),'units','norm','hori','center','vert','top');
    
end
set(hax,'xlim',[0 inf],'ylim',[0 inf],'ytick',[],'xtick',[],'visible','off');

if nax>1
    [mn, mx] = matchy(hax);
else
    mx = max([R{stimlist}]+[E{stimlist}]);
end

axes(hax(1));
text(1,1,num2str(round(mx)),'units','norm','hori','left','vert','middle');

%% description

desc = sprintf('%s Series %d Exp %d ',unit.id,unit.iseries, unit.iexp);
set(gcf,'menubar','none','numbertitle','off','name',desc)

