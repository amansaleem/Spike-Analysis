function u = ExptTraceShow( expt, gain, offset, show)
% Plots the tuning for traces in an expt
%
% u = ExptTraceShow loads the experiment from global PICK
%
% u = ExptTraceShow( expt ) loads the experiment in expt
%
% u = ExptTraceShow( expt, gain ) yada
%
% u = ExptTraceShow( expt, gain, offset) yada
%
% u = ExptTraceShow( expt, gain, offset, show) lets you specify whether to
% show the traces in a matrix of plots (show = 1) or not (show = 0,
% DEFAULT)
%
%USAGE: expt is specified or is PICK.expt from ExptPick
%       gain is gain of multiclamp software; i.e. 50 for 50mv/mv
%       offset is voltage outside of cell
% 2010-05 Matteo and Bilal
% 2010-06 MC added output of u.traces




% msgbox('ExptTraceShow is obsolete -- use ExptPick and click on Traces -> Inspect');





SetDefaultDirs;
global PICK

if nargin< 4
    show = 0;
end

if nargin<3
    offset = 0;
end

if nargin<2
    gain = 1;
end

if nargin < 1
    expt = [];
end

if isempty(expt)
    if ~isfield(PICK,'expt') || isempty(PICK.expt)
        PICK.expt = ExptLoad(PICK);
    end
    
    expt = PICK.expt;
end

if isempty(expt)
    fprintf('No data?\n');
    return
end

p = ProtocolLoad(expt);

% to look at the data

u.traces = cell(expt.nstim,expt.nrepeats);

if show ==1
figure; clf; ax = zeros(expt.nstim,expt.nrepeats);
else
end

for istim = 1:expt.nstim
    
    for irpt = 1:expt.nrepeats
        if show ==1
            ax(istim,irpt) = gridplot(expt.nstim, expt.nrepeats, istim, irpt);
        else
        end

        if ~isempty(expt.data{1}{istim, irpt})
%             vv = double(expt.data{1}{istim, irpt});   
            vv = double(expt.data{1}{istim, irpt}); % for filtered LFP on Axon secondary, AI6  BH 10.07.10
%             vv = (vv*3276.8)./1000; % takes care of unitspervolt in ExptLoad line 214; then converted to raw mV per V
            if vv(1,1)<=0
%                 vv = ((vv*gain) + (abs(nanmean(vv(1,1:15)*gain)))); %lfp
                vv = ((vv*gain) + offset); %vm
            else
%                 vv = ((vv*gain) - (abs(nanmean(vv(1,1:15)*gain)))); %lfp
                vv = ((vv*gain) + offset); %vm
            end
            u.traces{istim,irpt} = vv;
            tt = (1:length(vv))/(expt.samplerate);
            if show ==1
                plot( tt, vv ); hold on
            else
            end
        end
    end
    if show ==1
    sidetitle( sprintf('Stim %d',istim) );
    else 
    end
end

if show ==1
    set(ax,'ylim',[-inf inf]);
    matchy(ax);
    xlabel('Time (s)');
else
end

%set(ax(1:end-1),'xticklabel',[]);
%set(ax,'box','off');
