function rrr = UnitRateMatrix(u,ns)
% Bins the spikes of a unit into a 3-d matrix of firing rates
%
% rrr = UnitRateMatrix(u,ns)
% returns rrr with dimensions ns x nstim x nrepeats
%
% WARNING: it assumes all stimuli have same length!!
%
% 2012-07 MC

dur = median(u.stimdurs(:)); % assumes all stimuli have same length!!

bindur = dur / ns;

rrr = zeros( ns, u.nstim, u.nrepeats );
for istim = 1:u.nstim
    for irepeat = 1:u.nrepeats
        spiketimes = u.spiketimes{istim,irepeat};
        for ispike = 1:length(spiketimes)
            is = ceil(spiketimes(ispike)/bindur);
            if is<=ns
                rrr( is, istim, irepeat ) = rrr( is, istim, irepeat )+1;
            end
        end
    end
end
