function h = PopPlotRasters( uu, rptlist, stimlist, magicunit, magictime )
% PopPlotRasters plots the rasters of the responses of set of units
%
% PopPlotRasters( uu, rptlist, stimlist, magicunit, magictime )
%
% PopPlotRasters( uu, rptlist, stimlist ) lets you specify which repeats
% and which stimulus you want to see the population response for. Both
% default to 1 if left empty
%
% PopPlotRasters( uu, rptlist, stimlist, magicunit ) plots the unit you are
% interested in in red (the others are black)
%
% PopPlotRasters( uu, rptlist, stimlist, magicunit, magictime ) lets you
% specify a magic time t (in sec) where things were supposed to happen (for 
% example the time of electrical stimulation).
%
% EXAMPLE:
% uu = UnitLoad(DIRS.spikes,'M110306TS',4,2);
% PopPlotRasters(uu,[1:10],1,1,2);
%
% part of Spikes

% 2011-03 MC
% 2011-03 MS added the possibility of plotting either the response to one
%         stimulus over many repeats, or the response on one repeat to many
%         stimuli, and added magicunit
%
% BE CAREFUL: the function currently only supports either one stimulus, many
% repeats, or many repeats, one stimulus, not many repeats and many
% stimuli. 


if nargin<5
    magictime = NaN;
end

if nargin<4
    magicunit = NaN;
end

if isempty(stimlist)
    stimlist = 1;
end

if isempty(rptlist)
    rptlist = 1;
end

nu = length(uu);
nstim = uu(1).nstim;
nrpts = uu(1).nrepeats;

if any(rptlist<1) || any(rptlist>nrpts)
    error('Bad repeat number');
end

if any(stimlist<1) || any(stimlist>nstim)
    error('Bad stimulus number');
end

if length(rptlist) > 1
     triallist = rptlist;
     maxtrial = length(rptlist);
     titlename = ['Stimulus ' int2str(stimlist)]; 
elseif length(stimlist) > 1
     triallist = stimlist;
     maxtrial = length(stimlist);
     titlename = ['Repeat ' int2str(rptlist)]; 
else triallist = 1;
     maxtrial = 1;
     titlename = 'Repeat 1 and Stimulus 1'; 
end

ncols = ceil(maxtrial/20);
nrows = ceil(maxtrial/ncols);

h = figure;
rax = zeros(maxtrial,1); % raster axes
itrial = 1; 

for istim = stimlist
    
    for irep = rptlist
        
        rax(itrial)= subplot(nrows,ncols,itrial); % axes('position',get(hax(istim),'position'));
        text(0,0.5,num2str(triallist(itrial)),'hori','right','vert','middle','units','norm','fontsize',18);
        hold on;
        
        plot([magictime magictime],[0 nu], 'b');
        
        for iu = 1:nu
            xx = uu(iu).spiketimes{istim,irep};
            if size(xx,1) > 1  % ensure it is a row vector
                xx = xx(:)';
            end
            if ~isempty(xx)
                if iu == magicunit, col = [1 0 0]; else col = [0 0 0]; end
                p = plot( [1;1]*xx, [iu-1;iu]*ones(size(xx)), '-','color',col);
                set(p,'linewidth',1);
            end
            plot([1;1]*uu(iu).stimdurs(istim,irep),[iu-1;iu],'r-');
        end
        
        maxdur = max(max(uu(1).stimdurs(:,:)));
        text(maxdur,itrial+1,num2str(round(maxdur*100)/100),'hori','center','vert','top');
        
        itrial = itrial+1;
        
    end

end

set(rax,'visible','on','ytick',[],'xlim',[0 inf], 'ylim',[0,nu]);
set(rax,'xcolor','w','ycolor','w');
set(rax,'xtick',[]);

title(rax(1),[uu(1).animal ' Expt ' num2str(uu(1).iseries) '-' num2str(uu(1).iexp) titlename ' ALL UNITS']);

set(gcf,'papertype','A4','paperunits','centimeters','paperposition',[0.63452 0.63452 19.715 28.408]);
