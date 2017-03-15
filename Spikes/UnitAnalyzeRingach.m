function [figs, FourierTuning] = UnitAnalyzeRingach(units, samplesperdeg, deltatmax, runflags)
% UnitAnalyzeRingach analyzes responses to a Ringach experiment
%
% UnitAnalyzeRingach(units) lets you specify which unit(s) to analyze
%
% UnitAnalyzeRingach(units, samplesperdeg) lets you
% specify the number of samples per degree (DEFAULT: 2)
%
% UnitAnalyzeRingach(units, samplesperdeg, deltatmax)
% lets you specify the max deltat in ms (DEFAULT: 150)
%
% UnitAnalyzeRingach(units, samplesperdeg, deltatmax, runflags)
% lets you specify runflags, a structure with fields set to 0 or 1:
%     DoReceptiveFields -- to see the receptive fields [DEFAULT: 0]
%     ContrastReceptiveFields -- to map contrast rather than luminance [1]
%     DoMultipleTuningCurves -- to see tuning curves at different times [0]
%     DoPeakTuningCurves -- tuning curves at time of maximum response [0]
%     DoBlankCorrect -- fairly obsolete [0];
%     PlotFullReceptiveFields -- to plot rf for ori and sf [1]
%
% figs = AnalyzeRingach(...) returns handles to the figures that were
% generated
%
% Description of the algorithm:
% The algorithm sets time-zero when a spike occures, then it goes from dt-min
% milliseconds (usually 50 ms) AFTER the spike to dt-max milliseconds (~200 ms)
% BEFORE the spike, and computes, every dt milliseconds (~9 ms), the probability
% that a grating was presented with a given orientation. Thus, time line for the algorithm
% goes the opposite way than the time of the physical process (stimulus->spike).
% Accordingly, the computed kernel shows the probability of each stimulus
% orientation at different times, typically from -50 to +200 ms.
%
% WATCH OUT FOR TIMING: Since the algorithm works "backward" in time,
% the probability that a given orientation elicited a spike is
% spread over the time duration of the stimulus (typically 32 ms, ie. 4 frames
% per orientation). The best guess is that the middle frame is the one that evoked
% the spike, so the response should be dalyed by ~16 ms.
%
% Example:
% SetDefaultDirs
% units = UnitLoad(DIRS.spikes,'catz075',5,16,[],998);
% UnitAnalyzeRingach(units,1)
%
% See also GetRandomStimInfo, GetRespProb, RingachGetFourierTuning
%
% 2004-11 Matteo Carandini
% 2007-04 MC added blank correction
% 2007-06 MC added curve for phase, added output figs
% 2007-09 MC added Ori vs Time figure
% 2007-09 MC changed samplesperdeg to 2
% 2007-09 MC added "program flag" ContrastReceptiveFields
% 2007-10 LB freqs, oris, phases and contrasts are only determined from stims, not from blank stims (l. 300) to avoid crashes when blank contains a different number of e.g. phases
% 2008-12 AB fixed the ~16 ms bug
% 2009-02 ND added a structure that you enter that contains flags that determine analysis routines
% 2009-06 SK made it return a cell array of kernels (FourierTuning), as it only returned a single kernel (for the last unit in the unit list)
% 2009-09 AZ Turned AnalyzeRingach into UnitAnalyzeRingach, which take units as input instead of animal, iseries, iexp

global DIRS;
global PICK;
global DEMO; %#ok<NUSED>

figs = [];

if nargin < 4
    DoReceptiveFields = 0; % set this to 1 if you want to see the receptive fields
    ContrastReceptiveFields = 1; % set to 1 if you want to do maps of contrast rather than luminance
    DoMultipleTuningCurves = 0; % tuning curves at different times
    DoPeakTuningCurves = 0; % tuning curves at time delay of maximum response
    DoBlankCorrect = 0;
    PlotFullReceptiveFields = 0; % plot one kernel for each unit as a function of ori and sf
else
    DoReceptiveFields = runflags.DoReceptiveFields;
    ContrastReceptiveFields = runflags.ContrastReceptiveFields;
    DoMultipleTuningCurves = runflags.DoMultipleTuningCurves;
    DoPeakTuningCurves = runflags.DoPeakTuningCurves;
    DoBlankCorrect = runflags.DoBlankCorrect;
    PlotFullReceptiveFields = runflags.PlotFullReceptiveFields;
end

if ~exist('DEMO','var')
    error('Make sure there is a global called DEMO that is not empty');
end

if isempty('DEMO')
    error('Make sure the global DEMO is not empty');
end

if nargin < 3 || isempty(deltatmax)
    deltatmax = 150; % ms
end

if nargin < 2 || isempty(samplesperdeg)
    samplesperdeg = 5; % was 2
end

if nargin < 1
   if exist('PICK','var')
      animal  = PICK.animal;
      iseries = PICK.iseries;
      iexp    = PICK.iexp;
      units   = UnitLoad(DIRS.spikes,animal, iseries, iexp);
   end
else
   animal  = units(1).animal;
   iseries = units(1).iseries;
   iexp    = units(1).iexp;
end

protocol = ProtocolLoad(        animal, iseries, iexp);

% unitlist = [];
% if isempty(unitlist)
%     unitlist = 1:length(units);
% end

deltats =  ( -30:5:deltatmax )/1000;
%deltats = ( -10:5:deltatmax )/1000;
ndeltats = length(deltats);

myscreen = ScreenLogLoad(animal, iseries, iexp);

% if isempty(myscreen)% create a fictitious myscreen
%     myscreen.FrameRate = 124.8794; % Hz;
%     myscreen.Dist = 40;
% end

myscreen.PixelSize = 1/samplesperdeg;
myscreen.Xmax = 30 /myscreen.PixelSize; % sets height to 30 deg
myscreen.Ymax = 20 /myscreen.PixelSize; % sets width  to 20 deg

[ StimFrames, StimPars, ProbFrame, FrameSequence ] = GetRandomStimInfo(protocol, myscreen);

nu = length(units);

nrows = round(sqrt(nu));
ncols = ceil(nu/nrows);

FigOriTime = figure;
OriTimeAx = zeros(nu,1);
for iu = 1:nu
    OriTimeAx(iu) = subplot(nrows,ncols,iu); hold on
end

if DoReceptiveFields
    FigHartleyTransform= figure('position',get(FigOriTime,'position')-[0 420 0 0]);
    HartleyAx = zeros(nu,1);
    for iu = 1:nu
        % HartleyAx(iu) = gridplot(nrows,ncols,iu); hold on
        HartleyAx(iu) = subplot(nrows,ncols,iu); hold on % ND changed from gridplot to subplot
    end
end


for iunit = 1:nu

%     iunit = unitlist(iu);

    fprintf(1, 'Analyzing unit %d of %d\n', iunit, length(units));
    ProbRespGivenFrame = GetRespProb(protocol, units(iunit), myscreen, ProbFrame, FrameSequence, deltats);

    %---------------------------------------------------------------------
    %               BLANK CORRECTION
    %---------------------------------------------------------------------

    % this is a hack: all xfiles should return this field...
    if ~isfield(StimPars,'c')
        for istim = 1:protocol.nstim
            nframes = size( StimFrames{istim}, 3 );
            StimPars(istim).c = ones(nframes,1);
        end
    end

    if DoBlankCorrect

        % compute the probability of response to a blank stimulus
        MeanBlankResponse = zeros(1,length(deltats));
        nBlanks = 0;
        for istim = 1:protocol.nstim
            iiBlanks = (StimPars(istim).c == 0);
            nBlanks = nBlanks+ nnz(iiBlanks);
            MeanBlankResponse = MeanBlankResponse + sum( ProbRespGivenFrame{istim}(iiBlanks,:), 1 );
        end
        if nBlanks>0
            MeanBlankResponse = MeanBlankResponse / nBlanks;
        end

        % subtract it from the previous probabilities
        fprintf(1, 'Removing blank response\n');
        for istim = 1:protocol.nstim
            nframes = size( StimFrames{istim}, 3 );
            ProbRespGivenFrame{istim} = ProbRespGivenFrame{istim} - repmat( MeanBlankResponse, [nframes 1] );
        end

    end
    %-------------------------------------------------------------------
    %
    %                       the receptive fields
    %
    %-------------------------------------------------------------------

    if DoReceptiveFields

        % -------  get the mean StimFrame for each stimulus -------
        meanimage = cell(protocol.nstim,1);
        CorrectionImage = cell(protocol.nstim,1);
        for istim = 1:protocol.nstim
            fprintf(1,'Computing average image for stimulus %d\n', istim);
            nx = size(StimFrames{istim},1);
            ny = size(StimFrames{istim},2);
            meanimage{istim} = zeros(nx,ny,ndeltats);
            CorrectionImage{istim} = ones(nx,ny);
            nframes = size( StimFrames{istim}, 3 );
            for iframe = 1:nframes
                switch ContrastReceptiveFields
                    case 0
                        ThisFrame = StimFrames{istim}(:,:,iframe);
                    case 1
                        ThisFrame = abs(StimFrames{istim}(:,:,iframe));
                        CorrectionImage{istim} = CorrectionImage{istim}+ThisFrame;
                end
                % this is the key summation:
                for ideltat = 1:ndeltats
                    meanimage{istim}(:,:,ideltat) = ...
                        meanimage{istim}(:,:,ideltat) + ...
                        double(ProbRespGivenFrame{istim}(iframe,ideltat)*ThisFrame);
                end
            end
        end

        % ProbRespGivenFrame could be <0 because of blank correction

        % ---------------- average across stimuli ----------------
        StimPositions = cat(1,StimPars.position);
        minx = min(StimPositions(:,1));
        miny = min(StimPositions(:,2));
        maxx = max(StimPositions(:,3));
        maxy = max(StimPositions(:,4));

        nx = maxx - minx;
        ny = maxy - miny;

        grandsum = zeros(nx,ny,ndeltats);
        for istim = 1:protocol.nstim
            xx = ( StimPositions(istim,1):StimPositions(istim,3)-1 )-minx + 1;
            yy = ( StimPositions(istim,2):StimPositions(istim,4)-1 )-miny + 1;
            % hack added by MC 2007-09
            % xx = unique(ceil(xx/samplesperdeg));
            % yy = unique(ceil(yy/samplesperdeg));
            for ideltat = 1:ndeltats
                grandsum(xx,yy,ideltat) = grandsum(xx,yy,ideltat) + ...
                    meanimage{istim}(:,:,ideltat)./(eps+CorrectionImage{istim});
            end
        end
        grandsum = grandsum / max(max(max(abs(grandsum))));

        %------------ find the RF with maximal variance ------------
        vv = zeros(ndeltats,1);
        for ideltat = 1:ndeltats
            dd = grandsum(:,:,ideltat);
            vv(ideltat) = nanvar(dd(:));
        end
        [mx, myideltat] = max(vv);
        myideltat = myideltat(1); % just in case it finds two

        %----------------------------------------------------------------------
        %               plot all the RFs
        %----------------------------------------------------------------------

        figure(FigHartleyTransform);
        axes(HartleyAx(iunit));
        imagesc(squeeze(grandsum(:,:,myideltat)));
        axis equal; axis tight; axis ij; colormap bone
        set(gca,'xtick',[],'ytick',[]);
        title(sprintf('%s %2d ms',units(iunit).id,deltats(myideltat)*1000));

        figs(end+1) = figure; clf; ax = [];
        
        for ideltat = 1:ndeltats
            ax(ideltat) = subplot(ceil(ndeltats/8), 8, ideltat);
            imagesc( squeeze(grandsum(:,:,ideltat)) ); % , [ -1 1 ]
            xlabel(sprintf('%2d',deltats(ideltat)*1000));
        end
        set(ax,'clim',[min(grandsum(:)) max(grandsum(:))]);
        set(ax,'dataaspectratio',[1 1 1]);
        set(ax,'xtick',[], 'ytick',[]);
        colormap gray
        supertitle(sprintf('%s Expt %d-%d Unit %d-%d',...
            animal,iseries,iexp,units(iunit).ichan,units(iunit).icell))

        %----------------------------------------------------------------------
        %           plot the RF with maximal variance
        %----------------------------------------------------------------------

        %     xmindeg = ltpix2deg( minx - myscreen.Xmax/2, myscreen );
        %     ymindeg = ltpix2deg( miny - myscreen.Ymax/2, myscreen );
        %     xsizdeg = ltpix2deg( nx, myscreen );
        %     ysizdeg = ltpix2deg( ny, myscreen );
        %
        %    % compute center of mass
        %     [xmass,ymass] = meshgrid([1:nx],[1:ny]);
        %     zz = abs(grandsum(:,:,myideltat));
        %     xc = sum(sum(xmass.*zz))/sum(sum(zz));
        %     yc = sum(sum(ymass.*zz))/sum(sum(zz));
        %     xcdeg = ltpix2deg( minx+xc-myscreen.Xmax/2, myscreen );
        %     ycdeg = ltpix2deg( miny+yc-myscreen.Ymax/2, myscreen );
        %
        %     % compute peak
        %     [iy, ix] = find( zz == max(zz(:)) );
        %     iy = iy(1);
        %     ix = ix(1);
        %     xpdeg = ltpix2deg( minx+ix-myscreen.Xmax/2, myscreen );
        %     ypdeg = ltpix2deg( miny+iy-myscreen.Ymax/2, myscreen );
        %
        %     figs(end+1) = figure; clf
        %     imagesc(grandsum(:,:,myideltat))
        %     colormap gray
        %     set(gca,'dataaspectratio',[1 1 1])
        %     title(sprintf('%s Expt %d-%d Unit %d-%d. Delay = %2d ms',...
        %         animal,iseries,iexp,units(iunit).ichan,units(iunit).icell,deltats(myideltat)*1000));
        %     cmax = max(max(max(abs(grandsum))));
        %     caxis([-cmax,cmax]);
        %     grid on; hold on;
        %     plot(xc,yc,'bx');       % center of mass
        %     plot(ix,iy,'ro');       % peak
        %     set(gca,'xtick',linspace(0,nx,7),'xticklabel',round(10*linspace(xmindeg,xmindeg+xsizdeg,7))/10);
        %     set(gca,'ytick',linspace(0,ny,7),'yticklabel',round(10*linspace(ymindeg,ymindeg+ysizdeg,7))/10);
        %
        %     legend(...
        %         sprintf('(%2.2f,%2.2f) center of mass',xcdeg,ycdeg),...
        %         sprintf('(%2.2f,%2.2f) peak', xpdeg, ypdeg)        );
    end

    %----------------------------------------------------------------------
    %
    %                       The tuning curves
    %
    %----------------------------------------------------------------------

    % to avoid crashes when different number of phases etc. are specified
    % for the blank vs the true stims
    if protocol.nstim == protocol.blankstims
        allsfs  = [StimPars(1:end-1).sf ];
        alloris = [StimPars(1:end-1).ori];
        allphases = [StimPars(1:end-1).phase];
        allcs   = [StimPars(1:end-1).c];
    else
        allsfs = [StimPars.sf];
        alloris = [StimPars.ori];
        allphases = [StimPars.phase];
        allcs = [StimPars.c];
    end

    diffsf  = unique(allsfs(  allcs>0) ); diffsf ( isnan(diffsf ) ) = [];
    diffori = unique(alloris( allcs>0) ); diffori( isnan(diffori) ) = [];
    diffphase = unique(allphases( allcs>0 )); diffphase( isnan(diffphase) ) = [];

    ndiffsf  = length(diffsf);
    ndiffori = length(diffori);
    ndiffphase = length(diffphase);

    % Compute the tuning in Fourier space

    FourierTuning{iunit} = zeros(ndiffori,ndiffsf,ndeltats);

    for idiffsf = 1:ndiffsf
        for idiffori = 1:ndiffori
            count = 0;
            for istim = 1:protocol.nstim
                nframes = size( StimFrames{istim}, 3 );
                for iframe = 1:nframes
                    if  StimPars(istim).sf (iframe) == diffsf (idiffsf) && ...
                            StimPars(istim).ori(iframe) == diffori(idiffori)
                        count = count + 1;
                        for ideltat = 1:ndeltats
                            FourierTuning{iunit}(idiffori,idiffsf,ideltat) = FourierTuning{iunit}(idiffori,idiffsf,ideltat) + ProbRespGivenFrame{istim}(iframe,ideltat);
                        end % ideltat
                    end
                end % iframe
            end % istim
            FourierTuning{iunit}(idiffori,idiffsf,:) = FourierTuning{iunit}(idiffori,idiffsf,:)/count;
        end % idiffori
    end % idiffsf

    % Plot the RF in Fourier Space

    if PlotFullReceptiveFields

        if ndiffsf > 1 && ndiffori > 1
            figs(end+1) = figure;colormap('gray');
            ax = []; clf
            nax = ceil(sqrt(ndeltats));
            cmin = min(min(min(abs(FourierTuning{iunit}))));
            cmax = max(max(max(abs(FourierTuning{iunit}))));
            for it = 1:ndeltats
                ax(it) = subplot(nax,nax,it, 'align');
                imagesc(FourierTuning{iunit}(:,:,it));
                title(num2str(deltats(it)*1000));
                caxis([cmin,cmax]);
                set(gca,'xtick',[],'ytick',[]);
            end
            axes(ax(1));
            ylabel('Orientation');
            axes(ax(end));
            xlabel('Spatial frequency');
            supertitle(sprintf('%s Expt %d-%d Unit %d-%d',animal,iseries,iexp,units(iunit).ichan,units(iunit).icell));
        end

    end

    figure(FigOriTime);
    axes(OriTimeAx(iunit));
    OriVsTime = squeeze(mean(FourierTuning{iunit},2));
    imagesc(deltats,[diffori 180],OriVsTime([1:end 1],:));
    set(gca,'xlim',[-inf inf],'ylim',[-inf inf],'ytick',[0:45:180],'plotboxaspectratio',[2 1 1]);
    title(units(iunit).id);
    colorbar

    % Find the best deltat
    vv = zeros(ndeltats,1);
    for ideltat = 1:ndeltats
        dd = FourierTuning{iunit}(:,:,ideltat);
        isnotnan = setdiff(1:length(dd(:)),find(isnan(dd(:))));
        vv(ideltat) = var(dd(isnotnan));
    end
    [mx, myideltat] = max(vv);
    myideltat = myideltat(1); % just in case it finds two

    % Find the best ori, best sf and best deltat
    MeanOriTuning = mean(FourierTuning{iunit}(:,:,myideltat),2);
    [maxresp, myiori] = max(MeanOriTuning);
    tuning.sf = FourierTuning{iunit}(myiori,:,myideltat);
    [maxresp, myisf] = max(tuning.sf);
    tuning.ori = FourierTuning{iunit}(:,myisf,myideltat);

    tuning.isf      = myisf;
    tuning.iori     = myiori;
    tuning.diffsf   = diffsf;
    tuning.diffori  = diffori;

    %
    %       the tuning for phase at best time, best ori and best sf
    %

    PhaseTuning = zeros(ndiffphase,1);
    count = 0;
    for idiffphase = 1:ndiffphase
        for istim = 1:protocol.nstim
            nframes = size( StimFrames{istim}, 3 );
            for iframe = 1:nframes
                if  ...
                        StimPars(istim).sf (iframe) == diffsf (myisf) && ...
                        StimPars(istim).ori(iframe) == diffori(myiori) && ...
                        StimPars(istim).phase(iframe) == diffphase(idiffphase)
                    count = count + 1;
                    PhaseTuning(idiffphase) = PhaseTuning(idiffphase) + ProbRespGivenFrame{istim}(iframe,myideltat);
                end
            end % iframe
        end % istim

        PhaseTuning(idiffphase) = PhaseTuning(idiffphase)/count;
    end

    %----------------------------------------------------------------------
    %               plot the tuning curves at peak
    %----------------------------------------------------------------------

    if DoPeakTuningCurves
        figs(end+1) = figure; clf; ax = [];

        ax(1) = subplot(3,1,1);
        if length(diffori)>1
            plot( [ tuning.diffori tuning.diffori(1)+180 ], tuning.ori([1:end 1]), 'ko-', 'markerfacecolor', 'k' );
        end
        xlabel('Orientation (deg)');
        set(gca,'xlim',[0 180], 'xtick',[0:45:180]);

        ax(2) = subplot(3,1,2);
        if length(diffphase)>1
            plot( [diffphase 360], PhaseTuning([1:end 1]), 'ko-', 'markerfacecolor', 'k' );
        end
        set(gca,'xlim',[0 360],'xtick',0:45:360);
        xlabel('Phase (deg)');

        ax(3) = subplot(3,1,3);
        if length(diffsf)>1
            plot( tuning.diffsf, tuning.sf, 'ko-', 'markerfacecolor', 'k' );
        end

        xlabel('Spatial frequency (cpd)');

        set(ax,'ylim',[-inf inf],'plotboxaspectratio',[3 1 1 ]);

        set(gcf,'name',sprintf('%s Expt %d-%d Unit %d-%d. Delay = %2d ms',animal,iseries,iexp,units(iunit).ichan,units(iunit).icell,deltats(myideltat)*1000));
        supertitle(sprintf('%d-%d at %2d ms', units(iunit).ichan,units(iunit).icell,deltats(myideltat)*1000));
    end

    %---------------------------------------------------------------------
    %           tuning curves at different time intervals
    %---------------------------------------------------------------------

    if DoMultipleTuningCurves

        nax = ceil(sqrt(ndeltats));

        if ndiffsf > 1
            figs(end+1) = figure; clf; ax = [];
            for it = 1:ndeltats
                ax(it) = subplot(nax,nax,it, 'align');
                plot(log10(diffsf), mean(FourierTuning{iunit}(:,:,it),1),'ko-','markerfacecolor',[0 0 0]);
                text(0,0,num2str(deltats(it)*1000),'units','norm','hori','left','vert','bottom');
                axis tight
            end
            %         set(ax,'box','off','ylim',[0 inf]);
            %         matchy(ax,'bottom');
            set(ax,'box','off','ylim',[-inf inf]);
            matchy(ax);

            set(ax(1:end-1),'xtick',[]);
            set(ax(2:end),'ytick',[]);
            set(ax(end),'xtick',log10(lognums),'xticklabel',lognums);
            xlabel('Spatial frequency (cpd)');
            set(gcf,'name',sprintf('%s Expt %d-%d Unit %d-%d',animal,iseries,iexp,units(iunit).ichan,units(iunit).icell));
            title(units(iunit).icell);
        end

        if length(diffori)>1
            figs(end+1) = figure; clf; ax = [];
            for it = 1:ndeltats
                ax(it) = subplot(nax,nax,it, 'align'); hold on
                yy = mean(FourierTuning{iunit}(:,:,it),2);
                plot([diffori diffori(1)+180], yy([1:end 1]),'ko-','markerfacecolor',[0 0 0]);
                text(0,0,num2str(deltats(it)*1000),'units','norm','hori','left','vert','bottom');
                axis tight
            end
            %         set(ax,'box','off','ylim',[0 inf]);
            %         matchy(ax,'bottom');
            set(ax,'box','off','ylim',[-inf inf]);
            matchy(ax);

            set(ax(1:end-1),'xtick',[]);
            set(ax(2:end),'ytick',[]);
            set(ax(end),'xtick', 0:90:180 );
            xlabel('Orientation (deg)');
            set(gcf,'name',sprintf('%s Expt %d-%d Unit %d-%d',animal,iseries,iexp,units(iunit).ichan,units(iunit).icell));
            title(units(iunit).icell);
        end

        drawnow
    end % DoMultipleTuningCurves

end % loop on iunit

