function compareModelPred(model, es, t_low, t_high, t_norm)

% Program to see if the low / high situation is better predicted by the
% normal condition

% t = es.contrast==0.72 & es.traj~=0 & ~isnan(es.traj);

% First... is there a different in the mean raet between the conditions

figure('Position',[300 300 550 750]);

Y = es.spikeTrain;


[Y_norm EV_norm] = model.testMap(es.traj(t_norm), es.spikeTrain(t_norm,:));
[Y_high EV_high] = model.testMap(es.traj(t_high), es.spikeTrain(t_high,:));
[Y_low  EV_low]  = model.testMap(es.traj(t_low),  es.spikeTrain(t_low,:));

EV_norm(EV_norm<0) = 0;
EV_high(EV_high<0) = 0;
EV_low(EV_low<0) = 0;

k = (EV_high>0 | EV_low >0);
k_all = (EV_high>0 | EV_low >0) & EV_norm >0;

meanRate_norm  = nanmean(Y(t_norm,:),1);
meanRate_low  = nanmean(Y(t_low,:),1);
meanRate_high = nanmean(Y(t_high,:),1);

[distr, X] = hist(EV_high(k)-EV_low(k),30);
meanDiff = nanmean(EV_high(k) - EV_low(k));
semDiff = nansem(EV_high(k) - EV_low(k));

%% Plotting stuff
subplot(2,2,1)
plot(meanRate_low(k), meanRate_high(k),'k.');
display([num2str(nanmean(meanRate_high(k)-meanRate_low(k))) '+ -' num2str(nansem(meanRate_high(k)-meanRate_low(k)))]);

axis tight;
xlims = xlim;
ylims = ylim;
axis equal;
axis([-0.01 max(xlims(2), ylims(2)) -0.01 max(xlims(2), ylims(2))]);
line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
xlabel('Mean rate, low contrast');
ylabel('Mean rate, high contrast');


subplot(2,2,2)
plot(EV_low(k), EV_high(k),'r.')
axis tight;
xlims = xlim;
ylims = ylim;
axis equal;
axis([-0.001 max(xlims(2), ylims(2)) -0.001 max(xlims(2), ylims(2))]);
line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
xlabel('EV_{low \leftarrow normal}')
ylabel('EV_{high \leftarrow normal}')

subplot(2,1,2)

bar(X,distr,'EdgeColor','none','FaceColor','r')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
line([meanDiff meanDiff], ylim, 'color',[.5 0 0],'linewidth',2)
line([meanDiff-semDiff meanDiff-semDiff], ylim, 'color',[.5 0 0],'linewidth',2, 'linestyle','--')
line([meanDiff+semDiff meanDiff+semDiff], ylim, 'color',[.5 0 0],'linewidth',2, 'linestyle','--')

line([0 0], ylim, 'color','k','linewidth',1,'linestyle','--')
xlabel('EV_{high \leftarrow normal} - EV_{low \leftarrow normal}')
ylabel('Number of neurons')
