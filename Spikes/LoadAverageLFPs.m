function [lfp, fs] = LoadAverageLFPs(ExptTag,nShanks,nSites)
% Loads the LFP averaged across repeats for a Cerebus experiment
%
% [lfp, fs] = LoadAverageLFPs(ExptTag) returns a cell array
% {lfp{1},lfp{2},...} (one per stimulus). Each entry is a matrix nShanks x
% nSites x nt, where nt is the number of samples. ExptTag as usual must
% have fields animal, iseries, iexpt.
%
% LoadAverageLFPs(ExptTag,nShanks,nSites) lets you specify the number of
% shanks and number of sites. This is OBSOLETE.
%
% MC 2012-01-19

if nargin<3
    myProbe = Probe( ExptTag.animal, ExptTag.iseries );
    nShanks = myProbe.nShanks;
    nSites  = myProbe.nSites;
end


% this is buggy
% [~,exptfp] = ExptLoadCerebusTracesFPs(ExptTag,[],[],[],'LFPonly');

SetDefaultDirs;

AnimalDir  = fullfile(DIRS.data,ExptTag.animal);
SeriesDir  = fullfile(AnimalDir,num2str(ExptTag.iseries));
ExpDir     = fullfile(SeriesDir,num2str(ExptTag.iexp));
FileNameLFP   = fullfile(ExpDir, 'Expt_CerebusTraces_LFP.mat');

if ~exist(FileNameLFP,'file')
    disp('Cannot load it, sorry');
end

foo = load(FileNameLFP);
exptfp = foo.lfp;

p = ProtocolLoad(ExptTag);


%% average across repeats 

fs = exptfp.samplerate;
nrepeats = exptfp.nrepeats;

lfp = cell(p.nstim,1);

for istim = 1:p.nstim
    
    dur = min(exptfp.stimdurs(istim,:));
    nt = floor(dur*fs); % we will crop here -- maybe some are longer...
    
    lfp{istim} = zeros(nShanks,nSites,nt,nrepeats);
    for irepeat = 1:nrepeats
        for iShank = 1:nShanks
            channels = nSites*(iShank-1)+(1:nSites);
            lfp{istim}(iShank,:,:,irepeat) = exptfp.data{istim,irepeat}(channels,1:nt);
        end
    end
    lfp{istim} = mean(lfp{istim},4); % average across repeats
    
end
