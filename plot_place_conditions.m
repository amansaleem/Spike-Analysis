function plot_place_conditions(es, cell_list)

if nargin<2
    cell_list = 1:length(es.spikeIDs);
end

allCont = unique(es.contrast(~isnan(es.contrast)));
if length(allCont)==4
    cont = allCont(3);
elseif legnth(allCont)>1
    cont = allCont(end-1);
else
    cont = allCont;
end

t = es.traj>0 & es.smthBallSpd>7 & ~isnan(es.smthBallSpd)...
    & es.gain==1 & es.roomLength==1 & es.contrast==cont;

trl = es.traj>0 & es.smthBallSpd>7 & ~isnan(es.smthBallSpd)...
    & es.gain==1 & es.roomLength<1 & es.contrast==cont;
trh = es.traj>0 & es.smthBallSpd>7 & ~isnan(es.smthBallSpd)...
    & es.gain==1 & es.roomLength>1 & es.contrast==cont;

tgl = es.traj>0 & es.smthBallSpd>7 & ~isnan(es.smthBallSpd)...
    & es.gain<1 & es.roomLength==1 & es.contrast==cont;
tgh = es.traj>0 & es.smthBallSpd>7 & ~isnan(es.smthBallSpd)...
    & es.gain>1 & es.roomLength==1 & es.contrast==cont;

tcl = es.traj>0 & es.smthBallSpd>7 & ~isnan(es.smthBallSpd)...
    & es.gain==1 & es.roomLength==1 & es.contrast<cont & es.contrast>0;
tch = es.traj>0 & es.smthBallSpd>7 & ~isnan(es.smthBallSpd)...
    & es.gain==1 & es.roomLength==1 & es.contrast>cont;
tc0 = es.traj>0 & es.smthBallSpd>7 & ~isnan(es.smthBallSpd)...
    & es.contrast==0;

figure('Position', [232         193        1449         854]);

for icell = cell_list
    subplot(335)
    plot_subplot(es, t, icell)
    
    
    subplot(331)
    plot_subplot(es, tcl, icell)
    title('Contrast low')
    
    subplot(332)
    plot_subplot(es, tc0, icell)
    title('Contrast zero')
    
    subplot(333)
    plot_subplot(es, tch, icell)
    title('Contrast high')
    
    subplot(334)
    plot_subplot(es, trl, icell)
    title('Room length low')
    
    subplot(336)
    plot_subplot(es, trh, icell)
    title('Room length high');
    axis tight
    
    subplot(337)
    plot_subplot(es, tgl, icell)
    title('Gain low')
    
    subplot(339)
    plot_subplot(es, tgh, icell)
    title('Gain high')
    
    subplot(338)
    title(num2str(icell));
    axis off
    
    pause
end

    function plot_subplot(es, t, icell)
        
        trialIDs = unique(es.trialID(t));
        for itr = 1:length(trialIDs)
            es.trialID(es.trialID==trialIDs(itr)) = itr;
        end
        
        plot(es.trialID(t), es.traj(t),'color',[.5 .5 .5])
        hold on;
        plot(es.trialID(es.lick & t), es.traj(es.lick & t),'.','color',[.7 .7 .7])
        plot(es.trialID(es.spikeTrain(:,icell)>0 & t), es.traj(es.spikeTrain(:,icell)>0 & t),'ro')
        title(['Cell id: ' es.spikeIDs{icell}]);
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
        hold off;
        axis tight
    end
end