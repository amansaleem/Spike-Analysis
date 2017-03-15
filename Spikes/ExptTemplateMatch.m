function ExptTemplateMatch( expt, ichan )
% ExptTemplateMatch
%
% ExptTemplateMatch( expt, ichan )
%
% This is work in progress by Matteo, it does not result in the
% creation of any Units, just inspection by eye
%
% EXAMPLE:
% SetDefaultDirs
% expt = ExptLoad('catz082',7,14);
% expt = ExptLoad('catz025',2,3);
% expt = ExptLoad('catz027',6,1); % two cells
% expt = ExptLoad('catz035',1,2);
% ExptTemplateMatch( expt, ichan );
%
% 2009-02 MC created

global DIRS

chans = ChanLoad(DIRS.spikes,expt.animal,expt.iseries,expt.iexp);

chan = chans(ichan);

expt = ExptFilter(expt,chan);

interpfactor = 4;

%% Get the templates

% the prototypes
prots = chan.prots;
tt = chan.tt;

% interpolate them
Templates = resample(prots',interpfactor,1);
TemplateTimes = linspace(min(tt),max(tt),length(tt)*interpfactor);

% trim them
Templates( TemplateTimes<-0.25 | TemplateTimes > chan.mindur, : ) = [];
TemplateTimes ( TemplateTimes<-0.25 | TemplateTimes > chan.mindur ) = [];

[TemplateNSamples,nTemplates] = size(Templates);

TemplatePeakPos = zeros(nTemplates,1);
for iTemplate = 1:nTemplates
    [mx TemplatePeakPos(iTemplate)]=max(abs(Templates(:,iTemplate)));
end

% show them
figure; clf; ax = [];
for iTemplate = 1:nTemplates
    ax(iTemplate) = gridplot(nTemplates,1,iTemplate,1);
    ipeak = TemplatePeakPos(iTemplate);
    plot(TemplateTimes,Templates(:,iTemplate)); hold on
    plot(TemplateTimes(ipeak),Templates(ipeak,iTemplate),'ko');
    title(sprintf('Template for cell %d', chan.cellids(iTemplate) ));
    axis square
    axis tight
    box off
end
matchy(ax);

%% Template match

ErrMinValues = cell(nTemplates,expt.nstim,expt.nrepeats);

for is = 1:expt.nstim
    fprintf(1,'Stimulus %d of %d',is,expt.nstim);
    for ir = 1:expt.nrepeats
        fprintf(1,'.');
        vv = resample(double(expt.data{ichan}{is,ir}),interpfactor,1)';
        nt = length(vv);

        ErrValues = zeros(nTemplates,nt)*NaN;

        for iTemplate = 1:nTemplates
            % this can be sped up by using convolution...
            for it = 1:nt-TemplateNSamples
                ii = it+(1:TemplateNSamples);
                ErrValues(iTemplate,it) = norm(vv(ii)-Templates(:,iTemplate));
            end
            
            ii = 2:(nt-TemplateNSamples);
            ErrMinPosList = 1 + find( ...
                ErrValues(iTemplate,ii)<nanmedian(ErrValues(iTemplate,:)) & ...
                ErrValues(iTemplate,ii)<ErrValues(iTemplate,ii-1) & ...
                ErrValues(iTemplate,ii)<ErrValues(iTemplate,ii+1) );

            ErrMinValues{iTemplate,is,ir} = ErrValues(iTemplate,ErrMinPosList);
        end

    end
    fprintf(1,'\n');
end
    
%% Show distributions of error values

ErrThresh = NaN(nTemplates,1); % at first we don't know the threshold
for iTemplate = 1:nTemplates
    AllErrMinValues = [ErrMinValues{iTemplate,:,:}];
    figure; hist(AllErrMinValues,200);
    xlabel('Error amplitude');
    title(sprintf('Template for cell %d', chan.cellids(iTemplate) ));
end

%% Get thresholds from user

for iTemplate = 1:nTemplates
    ErrorThreshold = inputdlg(sprintf('Threshold for cell %d?', chan.cellids(iTemplate)) );
    ErrThresh(iTemplate) = str2double(ErrorThreshold{1});
end

%% Template Matching with assignment

SpikeSamples = cell(nTemplates,expt.nstim,expt.nrepeats);

vlims = [ min([expt.data{ichan}{:}]) max([expt.data{ichan}{:}]) ];
elims = [ min([ErrMinValues{:}]) max([ErrMinValues{:}]) ];

for is = 1:expt.nstim
    fprintf(1,'Stimulus %d of %d',is,expt.nstim);
    for ir = 1:expt.nrepeats
        fprintf(1,'.');
        vv = resample(double(expt.data{ichan}{is,ir}),interpfactor,1)';
        nt = length(vv);

        ErrValues = zeros(nTemplates,nt)*NaN;

        for iTemplate = 1:nTemplates
            % this can be sped up by using convolution...
            for it = 1:nt-TemplateNSamples
                ii = it+(1:TemplateNSamples);
                ErrValues(iTemplate,it) = norm(vv(ii)-Templates(:,iTemplate));
            end
            
            ii = 2:(nt-TemplateNSamples);
            ErrMinPosList = 1 + find( ...
                ErrValues(iTemplate,ii)<nanmedian(ErrValues(iTemplate,:)) & ...
                ErrValues(iTemplate,ii)<ErrValues(iTemplate,ii-1) & ...
                ErrValues(iTemplate,ii)<ErrValues(iTemplate,ii+1) );

            ErrMinValues{iTemplate,is,ir} = ErrValues(iTemplate,ErrMinPosList);
            SpikeSamples{iTemplate,is,ir} = ErrMinPosList(ErrMinValues{iTemplate,is,ir}<ErrThresh(iTemplate));
            
            if ~isnan(ErrThresh(iTemplate)) && (ir==1||ir==expt.nrepeats)
                ii = SpikeSamples{iTemplate,is,ir};
                figure; clf; ax = [];
                ax(1) = gridplot(2,1,1,1); plot(ErrValues(iTemplate,:)); hold on
                plot(ii,ErrValues(iTemplate,ii),'ko');
                ylabel('Error');
                set(gca,'ylim',elims);
                ax(2) = gridplot(2,1,2,1); plot(vv); hold on
                plot(ii+TemplatePeakPos(iTemplate),vv(ii+TemplatePeakPos(iTemplate)),'ko');
                set(gca,'ylim',vlims);
            end
        end

    end
    fprintf(1,'\n');
end

