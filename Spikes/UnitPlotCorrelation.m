function UnitPlotCorrelation( u1, u2, resolution )
% Cross correlation of two units
%
% UnitPlotCorrelation( u1, u2 )
% plots the cross-correlation of spike times of units u1 and u2 (obtained
% with UnitLoad)
%
% UnitPlotCorrelation( u1, u2, resolution ) lets you specify the resolution
% in Hz (DEFAULT: 1000)
%
% OBVIOUSLY this is not going to make sense (and will probably crash) if
% the two units are not from the same experiment.
%
% Example: 
% 
% choose an experiment with ExptLoad
% u1 = UnitLoad(DIRS.spikes, PICK.animal, PICK.iseries, PICK.iexp, 1001, 999);
% u2 = UnitLoad(DIRS.spikes, PICK.animal, PICK.iseries, PICK.iexp, 1013, 999);
% UnitPlotCorrelation( u1, u2 )
%
% 2010-03 MC and AZ started it
% 2012-07 MC removed need to load protocol, added resolution

% 
% u1 = UnitLoad(DIRS.spikes, PICK.animal, PICK.iseries, PICK.iexp, 1021, 999);
% u2 = UnitLoad(DIRS.spikes, PICK.animal, PICK.iseries, PICK.iexp, 1032, 999);

% p = ProtocolLoad(u1);

if nargin < 3
    resolution = 1000;
end

nlags = 50;

xcs = zeros(2*nlags+1,u1.nstim,u1.nrepeats);

for irepeat = 1:u1.nrepeats
    fprintf('.');
    for istim = 1:u1.nstim
        dur = u1.stimdurs(istim,irepeat);
        
        nt = ceil(dur*resolution);
        
        rr = zeros( nt, 2);
        
        rr( ceil(resolution*u1.spiketimes{istim,irepeat}), 1 ) = 1;
        rr( ceil(resolution*u2.spiketimes{istim,irepeat}), 2 ) = 1;
        
        [xcs(:,istim,irepeat), lags] = xcov( rr(:,1), rr(:,2), nlags, 'coeff' );
%         xc(isnan(xc)) = 0;
%         xc_mean = xc_mean + xc;
    end
end
% xc_mean = xc_mean / (u1.nstim*u1.nrepeats);
xcs(isnan(xcs)) = 0;
xc_mean = mean(mean(xcs,3),2);

fprintf('\n');

figure; 
plot( lags / resolution, xc_mean );
xlabel('Time (s)');
ylabel('Correlation coefficient');
title(sprintf('Correlation between units %s and %s in Series %d Experiment %d', ...
    u1.id,u2.id,u1.iseries,u1.iexp));

