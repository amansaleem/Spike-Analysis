function vv = ShowStimulusFrame( protocol, animal, iseries, iexp, iframe, istim, nographics )
% ShowStimulusFrame shows one frame of each stimulus in an experiment
%
% vv = ShowStimulusFrame( protocol, animal, iseries, iexp, iframe )
% vv{istim} has dimensions [ny,nx,3], use imshow(vv{istim})
%
% vv = ShowStimulusFrame( protocol, animal, iseries, iexp, iframe, istim )
% Allow to specify the stim (can be [] too). vv has dimensions [ny,nx,3]
%
% vv = ShowStimulusFrame( protocol, animal, iseries, iexp, iframe, istim, 'nographics' )
% vv = ShowStimulusFrame( protocol, animal, iseries, iexp, iframe, [], 'nographics' )
%
% 2004-12 Matteo Carandini
% 2007-08 AB added vv argout and istim/nographics argin

global DEMO
DEMO = 1;

ng = 1; if nargin<7, ng = 0; end

allstims = 0;
if nargin<6 | isempty(istim)
    allstims = 1;
end

myscreen = ScreenLogLoad(animal,iseries,iexp);
if isempty(myscreen)
    error('Cannot load a screen log for this experiment');
end

xfilefunc = protocol.xfile(1:end-2);
if ~ng; figure; end

if allstims
    for istim = 1:protocol.nstim
        stim = feval( xfilefunc, protocol.pars(:,istim), myscreen );
        vv{istim} = vsGetStimulusFrame(stim,myscreen, iframe);
        if ~ng; subplot(1,protocol.nstim,istim); imshow(vv{istim}); title(num2str(istim));end
    end
    supertitle(sprintf('Expt %s-%d-%d, frame %d (all stims)', animal, iseries, iexp, iframe));
else
    stim = feval( xfilefunc, protocol.pars(:,istim), myscreen );
    vv = vsGetStimulusFrame(stim,myscreen, iframe);
    if ~ng
        imshow(vv); title(num2str(istim))
        supertitle(sprintf('Expt %s-%d-%d, frame %d stim %d', animal, iseries, iexp, iframe, istim));
    end
end

