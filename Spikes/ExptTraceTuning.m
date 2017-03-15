function ExptTraceTuning( expt, gain, offset, TuningType )
% Plots the tuning for traces in an expt
% USAGE: ExptTraceTuning(expt,gain,offset,TuningType)
% expt can be specified, or is PICK.expt through ExptPick
% gain is gain of multiclamp software (must be entered)
% offset is voltage outside of cell
% TuningType is 1 for mean, 2 for Std, 3 for ... to be coded still
% 
%
% 2010-06 Matteo and Bilal
% 2010-06 Bilal fixed gain and offset

%msgbox('ExptTraceTuning is obsolete -- use ExptPick and click on Traces -> Inspect, then save and click on Tuning');

global PICK % may be useful

if nargin < 4
    TuningType = 1;
end

if nargin< 3
    offset = 0;
end
 
if nargin < 2
    gain = 10;
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

SetDefaultDirs;
p = ProtocolLoad(expt);

%% to look at the data

% for istim = 1:expt.nstim
%     
%     for irpt = 1:expt.nrepeats
%         %ax(istim,irpt) = gridplot(expt.nstim, expt.nrepeats, istim, irpt);
%         ax(istim) = gridplot(expt.nstim, 1, istim, 1);
%         if ~isempty(expt.data{1}{istim, irpt})
%             vv = double(expt.data{1}{istim, irpt});
%             vv = (vv*3276.8)/1000; % takes care of unitspervolt in ExptLoad line 214; then converted to raw mV per V
%             vv = ((vv*gain) + offset);
%             tt = (1:length(vv))/(expt.samplerate);
%             plot( tt, vv ); hold on
%             ylim([-70 -40]);
%         end
%     end
%     sidetitle( sprintf('Stim %d',istim) );
% end
% 
% %matchy(ax);
% xlabel('Time (s)');


%% extract tuning 

TuningOptions = {'mean','std'}; % let's add F1, F2 in the future...

mm = zeros(expt.nstim, expt.nrepeats);

for istim = 1:expt.nstim
    for irpt = 1:expt.nrepeats
       vv = double(expt.data{1}{istim, irpt});
            vv = (vv*3276.8)/1000; % takes care of unitspervolt in ExptLoad line 214; then converted to raw mV per V
            vv = ((vv*gain) + offset);
        if isempty(vv),
            mm(istim, irpt ) = NaN;
        end
        
        switch TuningOptions{TuningType}
            case 'mean'
                mm(istim, irpt ) = mean(vv);
            case 'std'
                mm(istim, irpt ) = std(vv);
            otherwise
                error('what????');
        end
    end
end

%% plot them

nonblanks = setdiff(1:expt.nstim, p.blankstims);

m = nanmean(mm,2);
e = nansem(mm')';

Vrest    = mean( m(p.blankstims) );
Vrest_se = norm( e(p.blankstims) );

figure; 
fillplot( [0 expt.nstim], [1 1]*(Vrest+Vrest_se), [1 1]*(Vrest-Vrest_se), 0.5*[1 1 1] ); hold on
errorbar( nonblanks, m(nonblanks), e(nonblanks), 'ko' ); hold on
plot( nonblanks, m(nonblanks), 'ko', 'markerfacec', 'k' );
ylabel(TuningOptions{TuningType});
xlabel('Stimulus number');












