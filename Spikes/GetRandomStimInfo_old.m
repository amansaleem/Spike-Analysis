function [StimFrames, StimPars, ProbFrame, FrameSequence ] = GetRandomStimInfo(prot, screeninfo, ResizeProportion, strOverwrite)
% GetRandomStimInfo information on a random stimulus
%
% [StimFrames, StimPars, ProbFrame, FrameSequence ] = GetRandomStimInfo(prot)
% obtains information about a random stimulus.
%
% [StimFrames, StimPars, ProbFrame, FrameSequence ] = GetRandomStimInfo(prot, screeninfo)
% lets you specify screen info (typically obtained with ScreenLogLoad, but with fewer pixels)
%
% The outputs are:
% - StimFrames, a tensor with each individual frame (an image)
% - StimPars, the parameters of each frame (e.g. orientation, spatial frequency,...)
% - ProbFrame, the probability of showing each frame
% - FrameSequence, the sequence of frames shown
%
% Example:
% SetDefaultDirs;
% p = ProtocolLoad('catz072',5,21);
% [StimFrames, StimPars, ProbFrame, FrameSequence ] = GetRandomStimInfo(p);
% % look at frame 4 of stimulus 2
% figure; imagesc(StimFrames{2}(:,:,4)); colormap bone; axis image
% % orientation of that frame
% StimPars(2).ori(4)
% % when was it shown:
% find(FrameSequence{2} == 4)
%
% See also AnalyzeRingach, GetRespProb
%
% 2004-11 Matteo Carandini
% 2007-09 RAF added a StimFrames resizing option (useful for saving memory)
% 2007-09 MC added some pre-allocations to stop matlab from complaining
% 2008-01 MC noted problems with computing StimFrames with new OpenGL stimuli
% 2008-03 LB added the computation of StimFrames for OpenGL stimuli
% 2008-03 LB added the possibility to load StimInfo from a *.mat file (to save time for computation for ogl stims)
% 2009-07 MC makes it load screeninfo if not provided

global DEMO %#ok<NUSED>


if ~exist('DEMO','var')
    error('Make sure there is a global called DEMO that is not empty');
end

if isempty('DEMO')
    error('Make sure the global DEMO is not empty');
end

if nargin < 4
    strOverwrite = '';
end

if nargin < 3
    ResizeProportion = 1;
end

if nargin < 2
    screeninfo = ScreenLogLoad(prot.animal, prot.iseries, prot.iexp);
    % hack things so that we take up less memory
    screeninfo.PixelSize = 1/5;
    screeninfo.Xmax = 30 /screeninfo.PixelSize; % sets height to 30 deg
    screeninfo.Ymax = 20 /screeninfo.PixelSize; % sets width  to 20 deg
end

%-----------------------------------------------------------------------
xfilefunction = prot.xfile(1:(end-2)); % drop the '.x'
if strcmp(xfilefunction,'visringlog')  %this was added by RAF 6/8/2005 to deal with the "size bug"
    xfilefunction = 'vringlog';
end

if exist(xfilefunction,'file')~=2
    error('The function %s is not in the path\n', xfilefunction);
end

istim = 1; % for preallocation
fprintf(1,'Make stimulus %d of %d\n',istim,prot.nstim);
if ~strcmp( strOverwrite, 'overwrite')
    % check if the stim info for the first stim already exists
end

try
    [OneStim, OneStimPars] = feval( xfilefunction,prot.pars(:,istim),screeninfo);
catch ME
    disp(ME.message);
    disp('Are you sure this is a random stimulus???');
    return
end
stim        = repmat( OneStim    , prot.nstim, 1);
StimPars    = repmat( OneStimPars, prot.nstim, 1);
for istim = 2:prot.nstim
    fprintf(1,'Make stimulus %d of %d\n',istim,prot.nstim);
    [ stim(istim), StimPars(istim) ] = feval( xfilefunction,prot.pars(:,istim),screeninfo);
    % will crash here if the xfile does not give two outputs
end

% free memory
clear OneStim;
clear OneStimPars;

% the following is simply not going to work for new xfiles based on OpenGL,
% where one specifies a rotation and the position of a window...
OpenGLStim = isfield(stim,'ori');
  
FrameSequence = cell(prot.nstim,1);
StimFrames    = cell(prot.nstim,1);
ProbFrame    = cell(prot.nstim,1);
fprintf('Computing Frame information');
for istim = 1:prot.nstim % these are the seeds
    fprintf('.');
    StimPars(istim).position = stim(istim).position;
    
    if ~OpenGLStim % the traditional stimuli
        FrameSequence{istim} = stim(istim).sequence.frames;
        nframes = length(stim(istim).frames{1}); % frames are unique stimuli shown 
        [nx, ny] = size(stim(istim).frames{1}{1});
        StimFrames{istim} = zeros(ceil(nx.*ResizeProportion), ceil(ny.*ResizeProportion), nframes); %Added by RF 09-2007
        for iframe = 1:nframes
            if ResizeProportion ~= 1    %If structure added by RF 09-2007
                StimFrames{istim}(:,:,iframe) = imresize(stim(istim).frames{1}{iframe},ResizeProportion);
            else
                StimFrames{istim}(:,:,iframe) = stim(istim).frames{1}{iframe};
            end
        end
    else % ogl stimuli
        % this is the frame sequence with:
        % - length: number of monitor frames presented (dur * framerate)
        % - contains: index to the unique stimulus shown
        % - stim is the seed
        FrameSequence{istim} = StimPars(istim).VirtualSequence; 
                
        nframes = length(StimPars(istim).iMovieFrame); % frames are unique stimuli shown 
        [nx, ny] = size(StimPars(istim).screen{1});
        
        StimFrames{istim} = zeros(ceil(nx.*ResizeProportion), ceil(ny.*ResizeProportion), nframes, 'uint8');
        for iframe = 1:nframes
            if ResizeProportion ~= 1  
                StimFrames{istim}(:,:,iframe) = imresize(StimPars(istim).screen{iframe},ResizeProportion);
            else
                StimFrames{istim}(:,:,iframe) = StimPars(istim).screen{iframe};
            end
            StimPars(istim).screen{iframe} = []; % free some memory
        end
    end % ogl stimuli

    StimFrames{istim} = (single(StimFrames{istim}) - 3)/252; % now it is bet  0 and 1 % convert to single to save memory
    StimFrames{istim} = (single(StimFrames{istim}) - 0.5)*2; % now it is bet -1 and 1
    
    nseq = length(FrameSequence{istim});
    ProbFrame{istim} = zeros(1,nframes);
    for iframe = 1:nframes
        ProbFrame{istim}(iframe) = sum(FrameSequence{istim} == iframe)/nseq;
    end
end
fprintf('\n');

% this is a hack: all xfiles should return this field...
if ~isfield(StimPars,'c')
    warning('xfile %s does not return Stimpars.c -- setting it to 1 arbitrarily', prot.xfile);
    for istim = 1:prot.nstim
        nframes = size( StimFrames{istim}, 3 );
        StimPars(istim).c = ones(nframes,1);
    end
end

if ~isfield(StimPars,'sf')
    warning('xfile %s does not return Stimpars.sf -- setting it to 1 arbitrarily', prot.xfile);
    for istim = 1:prot.nstim
        nframes = size( StimFrames{istim}, 3 );
        StimPars(istim).sf = ones(1,nframes);
    end
end


