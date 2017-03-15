function [S,Savg,Sstd,rr] = UnitGetSparseness(u,p,myscreen,stim,icells)
% UNITGETSPARSENESS: Analyze natural movie stimuli for sparseness
% x.file used must have been vismovie
% S    = cell of sparseness indices, for each stimulus and repeat
% Savg = S averaged   across repeats, not including 0 spike runs
% Sstd = std dev of S across repeats, not including 0 spike runs
% rr   = number of spikes, binned by frame
% 
% Examples: 
% 
% SetDefaultDirs;
% u = UnitLoad(DIRS.spikes,'CATZ077',4,28,[],1);
% [S,Savg,Sstd,rr] = UnitGetSparseness(u);
% 
% - or -
% 
% SetDefaultDirs;
% u = UnitLoad(DIRS.spikes,'CATZ077',4,28,[],1);
% p = ProtocolLoad(        u(1).animal,u(1).iseries,u(1).iexp);
% myscreen = ScreenLogLoad(u(1).animal,u(1).iseries,u(1).iexp);
% stim = makestims(p,myscreen,1);
% [S,Savg,Sstd,rr] = UnitGetSparseness(u,p,myscreen,stim{1}{1});
% 
% 
% Created AZ 2009-10-10
% TODO: get rid of icells, uIX variables?
% TODO: what to do about 0 spike situations?

%% Load expt infos
if ~exist(       'p','var')
   p = ProtocolLoad(u(1).animal,u(1).iseries,u(1).iexp);
end
if ~strcmp(p.xfile,'vismovie.x')
   error('Protocol used for this expt was not vismovie.x.');
end

if ~exist('myscreen','var')
   myscreen = ScreenLogLoad(p.animal,p.iseries,p.iexp);
end
if ~exist(    'stim','var')
   stim = makestims(p,myscreen,1);
   stim = stim{1}{1};
end
if exist(   'icells','var')
   uIX = find(ismember([u.icell]',icells))';
else
   uIX = find([u.icell]')';
end
[nstims nreps] = size(u(uIX(1)).spiketimes);

numframes    =  max(stim.sequence.frames);
numframereps = mean(diff(find(~isnan(stim.sequence.frames))));
framedur     = numframereps / myscreen.RealFrameRate;

bins = 0:framedur:numframes*framedur;
clear myscreen stim numframereps framedur

%% Initialize, calculate sparseness
rr   = repmat({zeros(nstims,nreps)},numel(uIX),1);
S    = rr;
Savg = zeros(numel(uIX),nstims);
Sstd = zeros(numel(uIX),nstims);

for cellIX = 1:numel(uIX)
   % Per unit
   % u(16).spiketimes{1,1};

   % Calulcate num spikes in each time bin FOR ALL RUNS, & Sparseness
   rr{cellIX} = cellfun(@(x) histc(x,bins), u(uIX(cellIX)).spiketimes,'UniformOutput',false);
   A  = cellfun(@(x) (sum(x/numframes)^2) / sum((x.^2)/numframes),rr{cellIX});
   S{cellIX} = (1 - A) / (1 - 1/numframes);

   %% TODO: WHAT TO DO WITH 0 spike CASE?
%    expt.Savg(cellIX,:) = mean(expt.S{cellIX},2);
%    expt.Sstd(cellIX,:) =  std(expt.S{cellIX},[],2);
   for istim = 1:nstims
      Savg(cellIX,istim) = mean(S{cellIX}(istim,(~isnan(S{cellIX}(istim,:)))));
      Sstd(cellIX,istim) =  std(S{cellIX}(istim,(~isnan(S{cellIX}(istim,:)))));
   end
   % figure; bar(bins,rr{1,1},'histc')
end
clear istim cellIX

end