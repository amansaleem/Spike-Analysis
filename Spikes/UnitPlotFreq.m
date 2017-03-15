function fh = UnitPlotFreq( unit, stim, f, stimNum)
% UnitPlotFreq plots power spectrum of responses in a Unit data structure
%
% fh = UnitPlotFreq( unit )
%   returns the figure handle
%
% fh = UnitPlotFreq( unit, stim) allows to specify the stimuli for which to plot
% the power spectrum, e.g. [4 7 9; 5 8 10]
%
% fh = UnitPlotFreq( unit, stim, f) allows to specify the harmonic f of the
% stimulus to plot, e.g. f=2 would be the second harmonic
%
% fh = UnitPlotFreq( unit, stim, f, stimNum) allows to specify the component of
% the stimulus for which the frequency should be marked by a vertical bar
% (default: 1). Set to 2 if you want to see the freq of the 2nd orientation
% in a plaid.
%
% part of Spikes
%
% 12/00 MC
% 03-01-2002 VM, now uses the true tfs of the stimuli when possible
% 2003-03 VM made it fit for traces
% 2003-03 VM changed temporal resolution from 5 to 1 ms
% 2008-06 LB added input argument stim
% 2008-07 LB added input argument f, stimNum
% 2008-07 LB replaced search for stim freqs with function UnitGetStimFreqs
%   improved the plot
% 2008-07 LB added the output argument fh
%
% unit = UnitLoad( DIRS.spikes, 'CATZ008', 5, 2, 1, 2);

if nargin < 4
    stimNum = 1;
end
if nargin < 3
    f = 1;
end
if nargin < 2
    stim = [];
end

% Decide what datatype to use
datatype = unit.datatype; % can be 'spiketimes' or 'traces'

%data = getfield(unit,datatype);

if isempty(stim)
    nstim = size(unit.(datatype),1);
    stim = 1:nstim;
end
nrows = size(stim, 1);
ncols = size(stim, 2);

%nrpts	= size(unit.(datatype),2);

switch datatype
case 'spiketimes'
   [R, E, resolution, dt] = UnitGetRates( unit, 0.00025, 'gauss');
case 'traces'
   [R, E, resolution, dt] = UnitGetRates(unit);
end

%% figure

fh = figure; clf
ax = zeros(nrows,ncols);
for irow = 1 : nrows
    for icol = 1 : ncols
   
        ax(irow,icol)= gridplot(nrows,ncols,irow,icol); % axes('position',get(hax(istim),'position'));
        text(0,0.5,num2str(stim(irow,icol)),'hori','right','vert','middle','units','norm','fontsize',18);
        hold on;
   
        nsamples = length(R{stim(irow,icol)});
        duration = nsamples*dt;
   
        absrft = abs(fft(R{stim(irow,icol)}))/(nsamples/2);
        % we care for freqs up to 100 Hz,
        ff = freq(nsamples,duration);
        ii = find(ff>=0.5 & ff<=50);
        plot( ff(ii), absrft(ii));
   
        % OLD: indmax = findmax(absrft(ii).^2);
        [mymax, indmax] = max(absrft(ii).^2);
        if length(indmax)==1 && absrft(ii(indmax))> 5*mean(absrft(ii))
            text( ff(ii(indmax)), absrft(ii(indmax)), ...
                num2str(ff(ii(indmax)),4), 'hori', 'center','vert', 'bottom');
        end
        zoom on;
    end
end
set(ax,'visible','on','ylim',[0 inf],'xlim',[0.5 50],'xscale','log','yscale','linear');
set(ax,'xtick',[1 2 5 10 20 50],'xticklabel',[],'box','off');
set(ax(end),'xticklabel',[1 2 5 10 20 50]);
set(ax,'yaxislocation','right');

[mn, mx] = matchy(ax,'bottom');

xlabel(ax(end),'Frequency (Hz)');

desc = sprintf('Unit %s Exp %d - %d ',unit.id, unit.iseries, unit.iexp);
set(gcf,'numbertitle','off','name',desc)

%% Add vertical bars at the freqs that are in the stimulus

protocol = ProtocolLoad(unit.animal, unit.iseries, unit.iexp);
StimFreq = UnitGetStimFreqs(protocol, unit, stimNum);

% protocol = ProtocolLoad( unit.animal, unit.iseries, unit.iexp, 'loadscreen' );
% StimFreq = [];
% if isfield(protocol,'estfreqs') && ~isempty(protocol.estfreqs);
%    StimFreq = protocol.estfreqs;
% elseif isfield(protocol,'pfilefreqs')
%    StimFreq = protocol.pfilefreqs;
%    disp('WARNING: the stimulus temporal frequencies shown in the plot are NOT exact');
% elseif isfield(protocol,'ncycles') ~any(isnan(protocol.ncycles))
%    durs_data = mean(unit.stimdurs,2);
%    StimFreq = protocol.ncycles(:)./durs_data(:);
%    disp('WARNING: the stimulus temporal frequencies shown in the plot are NOT exact');
% end

if ~isempty(StimFreq)
   StimFreq = StimFreq(:);
   plotStim = setdiff(stim(:),protocol.blankstims);
   for istim = 1 : length(plotStim)
      myax = (plotStim(istim) == stim);
      axes(ax(myax));
      plot(f*StimFreq(plotStim(istim))*[1 1],[ 0 mx ], 'r:');
      text(f*StimFreq(plotStim(istim)), mx+5, sprintf('%2.3f', f*StimFreq(plotStim(istim))), 'Color', 'r');
   end
end
