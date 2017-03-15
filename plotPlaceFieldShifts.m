model_hRL = con.FigObjects.figGUIobj2.maps1d.highRoomLength;
model_n = con.FigObjects.figGUIobj2.maps1d.normal;
model_lg = con.FigObjects.figGUIobj2.maps1d.lowGain;
model_lRL = con.FigObjects.figGUIobj2.maps1d.lowRoomLength;
model_hg = con.FigObjects.figGUIobj2.maps1d.highGain;
model_lc = con.FigObjects.figGUIobj2.maps1d.lowContrast;
model_hc = con.FigObjects.figGUIobj2.maps1d.highContrast;


[meanModel_n, bins_n, cellList_n, maxPos] = getPopNormVector(model_n);
% Contrast
[meanModel_hc, bins_hc, cellList_hc, maxPos_hc] = getPopNormVector(model_hc, cellList_n,1);
[meanModel_lc, bins_lc, cellList_lc, maxPos_lc] = getPopNormVector(model_lc, cellList_n,1);
% Gain
% [meanModel_hg, bins_hg, cellList_hg, maxPos_hg] = getPopNormVector(model_hg, cellList_n,1);
% [meanModel_lg, bins_lg, cellList_lg, maxPos_lg] = getPopNormVector(model_lg, cellList_n,1);
% RoomLength
[meanModel_hRL, bins_hRL, cellList_hRL, maxPos_hRL] = getPopNormVector(model_hRL, cellList_n,1);
[meanModel_lRL, bins_lRL, cellList_lRL, maxPos_lRL] = getPopNormVector(model_lRL, cellList_n,1);

% maxPos_n = maxPos;

figure('Position',[680 430 850 700]);
% Contrast column
subplot(331)
imagesc(bins_lc, 1:size(meanModel_lc,1), meanModel_lc); axis xy;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
hold on; plot(bins_n(maxPos(cellList_lc)),1:length(cellList_lc),'r','linewidth',2); hold off;
title('Low Contrast')
subplot(334)
imagesc(bins_n, 1:size(meanModel_n,1), meanModel_n); axis xy;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
hold on; plot(bins_n(maxPos(cellList_n)),1:length(cellList_n),'r','linewidth',2); hold off;
title('Normal')
subplot(337)
imagesc(bins_hc, 1:size(meanModel_hc,1), meanModel_hc); axis xy;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
hold on; plot(bins_n(maxPos(cellList_hc)),1:length(cellList_hc),'r','linewidth',2); hold off;
title('High Contrast')
xlabel(' Position')
ylabel('Cell')

% % Gain column
% %figure;
% subplot(332)
% imagesc(bins_lg, 1:size(meanModel_lg,1), meanModel_lg); axis xy;
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% hold on; plot(bins_n(maxPos(cellList_lg)),1:length(cellList_lg),'r','linewidth',2); hold off;
% hold on; plot(bins_hRL(maxPos(cellList_lg)),1:length(cellList_lg),'g','linewidth',2); hold off;
% title('Low Gain')
% subplot(335)
% imagesc(bins_n, 1:size(meanModel_n,1), meanModel_n); axis xy;
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% hold on; plot(bins_n(maxPos(cellList_n)),1:length(cellList_n),'r','linewidth',2); hold off;
% title('Normal')
% subplot(338)
% imagesc(bins_hg, 1:size(meanModel_hg,1), meanModel_hg); axis xy;
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% hold on; plot(bins_n(maxPos(cellList_hg)),1:length(cellList_hg),'r','linewidth',2); hold off;
% hold on; plot(bins_lRL(maxPos(cellList_hg)),1:length(cellList_hg),'g','linewidth',2); hold off;
% title('High Gain')
% xlabel(' Position')
% ylabel('Cell')

% Room Length column
%figure;
subplot(333)
imagesc(bins_lRL, 1:size(meanModel_lRL,1), meanModel_lRL); axis xy;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca,'xlim',[0 max(bins_hRL)]);
hold on; plot(bins_n(maxPos(cellList_lRL)),1:length(cellList_lRL),'r','linewidth',2); hold off;
hold on; plot(bins_lRL(maxPos(cellList_hRL)),1:length(cellList_hRL),'g','linewidth',2); hold off;
title('Low RL')
subplot(336)
imagesc(bins_n, 1:size(meanModel_n,1), meanModel_n); axis xy;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
hold on; plot(bins_n(maxPos(cellList_n)),1:length(cellList_n),'r','linewidth',2); hold off;
set(gca,'xlim',[0 max(bins_hRL)]);
title('Normal')
subplot(339)
imagesc(bins_hRL, 1:size(meanModel_hRL,1), meanModel_hRL); axis xy;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
hold on; plot(bins_n(maxPos(cellList_hRL)),1:length(cellList_hRL),'r','linewidth',2); hold off;
hold on; plot(bins_hRL(maxPos(cellList_hRL)),1:length(cellList_hRL),'g','linewidth',2); hold off;
set(gca,'xlim',[0 max(bins_hRL)]);
title('High RL')
xlabel(' Position')
ylabel('Cell')

colormap(gray)