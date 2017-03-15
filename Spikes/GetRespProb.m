function ProbRespGivenFrame = GetRespProb(protocol, unit, myscreen, ProbFrame, FrameSequence, deltats)
% GetRespProb computes the conditional probability of a response in a noise experiment
%
% ProbRespGivenFrame = GetRespProb(protocol, unit, myscreen, ProbFrame, FrameSequence, deltats)
%
% See also AnalyzeRingach, GetRandomStimInfo
%
% 2004-11 Matteo Carandini
% 2007-06 MC added warning to deal with missing repeats
% 2007-09 MC set it to use RealFrameRate, and cleaned up
% 2008-12 AB fixed the 16ms bug

if isfield( unit, 'spiketimes' )
    datatype = 'delta';
    % unit is a structure of type "unit"
else
    datatype = 'continuous';
    % unit is a structure of type "expt" (the result of a model)
end

if size(unit.spiketimes, 2) < unit.nrepeats
    warning('Not all repeats have been discriminated');
    unit.nrepeats = size(unit.spiketimes, 2);
end

%%

ProbFrameAndResp    = cell(protocol.nstim,1);
ProbRespGivenFrame  = cell(protocol.nstim,1);

ndeltats = length(deltats);

%----- Added to fix the 16 ms problem - AB Dec 14, 2008-------
nfrmsProtocol = protocol.pars(strmatch('nfr',protocol.parnames),1);
nfrmsNoise    = min(diff(find(diff(FrameSequence{1}))));
if nfrmsProtocol~=nfrmsNoise
    error('Cannot determine the number of frames per orientation!');
else
    delay = (nfrmsProtocol/myscreen.RealFrameRate)/2;
end
%-------------------------------------------------------------

for istim =1:protocol.nstim  
    
    fprintf(1,'Analyzing stimulus %d of %d\n',istim,protocol.nstim);
    nseq = length(FrameSequence{istim});
    nframes = length(ProbFrame{istim});
    ProbFrameAndResp{istim} = zeros(nframes,ndeltats);
   
    for ideltat = 1:ndeltats          
        deltat = deltats(ideltat);
        
        for irepeat = 1:unit.nrepeats
            
            switch datatype
                
                case 'delta'
                    
                    % the frames presented deltat seconds before spikes
                    % indexlist = floor( (unit.spiketimes{istim,irepeat} - deltat) * myscreen.FrameRate );
                    % AB added 'delay' to fix 16ms problem (Dec 14, 2008)
                    indexlist = floor( (unit.spiketimes{istim,irepeat} - deltat + delay) * myscreen.RealFrameRate );
                    
                    indexlist(indexlist<=0 | indexlist>nseq) = [];
                    if length(indexlist)>1
                        framelist = FrameSequence{istim}( indexlist );
                        frameprobs = hist(framelist, 1:nframes )/nseq;
                        % averages across repeats
                        ProbFrameAndResp{istim}(:,ideltat) = ProbFrameAndResp{istim}(:,ideltat) + frameprobs'/unit.nrepeats;            
                    end
                    
                case 'continuous'
                    
                    resp = unit.data{1}{istim,irepeat};
                    nsamples = length(resp);
                    
                    deltasamples = round( deltat*unit.samplerate + 0.5);
                    
                    % new code (by Matteo)
                    ShiftedResp = resp( mod( (1:nsamples) + deltasamples - 1, nsamples ) + 1 );
                    ShiftedResp(nsamples-deltasamples:nsamples) = NaN;
                    ShiftedResp(1:(-deltasamples)) = NaN;
                    
                    % old code (by Valerio)
%                     respIndex = max(deltasamples,1):nsamples;
%                     if deltasamples < 1
%                         ShiftedResp(1:1-deltasamples) = 0;
%                         ShiftedResp(2-deltasamples:nsamples) = resp(respIndex(1:length(respIndex)+deltasamples-1));
%                     else
%                         ShiftedResp(1:length(respIndex)) = resp(respIndex);
%                         ShiftedResp(length(respIndex)+1:nsamples) = 0;
%                     end            
                    
                    for iframe = 1:nframes
                        % tt = find(FrameSequence{istim} == iframe)/myscreen.FrameRate;
                        tt = find(FrameSequence{istim} == iframe)/myscreen.RealFrameRate;
                        ii = round(tt*unit.samplerate);
                        ProbFrameAndResp{istim}(iframe,ideltat) = nanmean(ShiftedResp(ii));
                    end
            end % switch
            
        end % loop on irepeat
        
        weights = ProbFrame{istim}(:);
        weights( weights==0 ) = NaN;
        ProbRespGivenFrame{istim}(:,ideltat) = ProbFrameAndResp{istim}(:,ideltat)./weights;
        
        foo = ProbRespGivenFrame{istim}(:,ideltat);
        foo( ~isfinite(foo) ) = 0;
        ProbRespGivenFrame{istim}(:,ideltat) = foo;
        
    end % loop on ideltat
    
end % loop on istim

