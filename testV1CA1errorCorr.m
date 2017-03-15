% % load Posterior_all_150226_trainedCorrect_gaussSmth_250_noInter
% load Posterior_all_150327_trainedCV_250_noInter
% CA = Posterior_all;
% CA = CA([1:6 8]);
% % % 
% % load Posterior_all_151014_V1_trainedCV_125
% load Posterior_all_V1_150601_trainedCV_250.mat
% % load Posterior_all_151013_V1_trainedCorrect_250_quick
% V1 = Posterior_all;
clear delay_* V1CA1* peVC pVC Posterior*

thres = 10;

start = 0;
stop  = 50;
for n = 1:7%[1:3 5:7]
    % Low Contrast
    k_l = CA(n).X_low_orig>start & CA(n).X_low_orig<=stop;
%     outcome = CA(n).data.outcome(CA(n).t_low)~=2;
%     k_l = outcome(1:length(k_l)) & k_l; 
    C_err_l = CA(n).error.low(k_l)';
    V_err_l = V1(n).error.low(k_l)';
    sub_l = 0;%25-CA(n).X_low_orig(k_l);
    [V1CA1corr(n,1), pVC(n,1)] = corr(C_err_l-sub_l, V_err_l-sub_l);
    [delay_corr_l(n,:), lags] = xcov(C_err_l-sub_l, V_err_l-sub_l, 30);
    
%     [V1CA1ecorr(n,1), peVC(n,1)] = corr(sign(V_err_l).*(abs(V_err_l)>thres), sign(C_err_l).*(abs(C_err_l)>thres));
%     [delay_ecorr_l(n,:), elags] = xcov(sign(V_err_l).*(abs(V_err_l)>thres), sign(C_err_l).*(abs(C_err_l)>thres), 30);
    [V1CA1ecorr(n,1), peVC(n,1)] = corr((abs(V_err_l)>thres), (abs(C_err_l)>thres));
    [delay_ecorr_l(n,:), elags] = xcov((abs(V_err_l)>thres), (abs(C_err_l)>thres), 30);
%     [V1CA1ecorr(n,1), peVC(n,1)] = corr((abs(V_err_l)>thres), (abs(C_err_l(randperm(length(C_err_l))))>thres));
%     [delay_ecorr_l(n,:), elags] = xcov((abs(V_err_l)>thres), (abs(C_err_l)>thres), 30);
    
    % Norm Contrast
    k_n = CA(n).X_norm>start & CA(n).X_norm<=stop;
%     outcome = CA(n).data.outcome(CA(n).t_norm)~=2;
%     k_n = outcome(1:length(k_n)) & k_n; 
    C_err_n = CA(n).error.norm(k_n)';
    V_err_n = V1(n).error.norm(k_n)';
    sub_n = 0;%25-CA(n).X_norm(k_n);
    [V1CA1corr(n,2), pVC(n,2)] = corr(C_err_n-sub_n, V_err_n-sub_n);
    [delay_corr_n(n,:), lags] = xcov(C_err_n-sub_n, V_err_n-sub_n, 30);
    
%     [V1CA1ecorr(n,2), peVC(n,2)] = corr(sign(V_err_n).*(abs(V_err_n)>thres), sign(C_err_n).*(abs(C_err_n)>thres));
%     [delay_ecorr_n(n,:), elags] = xcov(sign(V_err_n).*(abs(V_err_n)>thres), sign(C_err_n).*(abs(C_err_n)>thres), 30);
    [V1CA1ecorr(n,2), peVC(n,2)] = corr((abs(V_err_n)>thres), (abs(C_err_n)>thres));
    [delay_ecorr_n(n,:), elags] = xcov((abs(V_err_n)>thres), (abs(C_err_n)>thres), 30);
%     [V1CA1ecorr(n,2), peVC(n,2)] = corr((abs(V_err_n)>thres), (abs(C_err_n(randperm(length(C_err_n))))>thres));
%     [delay_ecorr_n(n,:), elags] = xcov((abs(V_err_n)>thres), (abs(C_err_n)>thres), 30);
    
    % High Contrast
    k_h = CA(n).X_high_orig>start & CA(n).X_high_orig<=stop;
%     outcome = CA(n).data.outcome(CA(n).t_high)~=2;
%     k_h = outcome(1:length(k_h)) & k_h; 
    C_err_h = CA(n).error.high(k_h)';
    V_err_h = V1(n).error.high(k_h)';
    sub_h = 0;%25-CA(n).X_high_orig(k_h);
    [V1CA1corr(n,3), pVC(n,3)] = corr(C_err_h-sub_h, V_err_h-sub_h);
    [delay_corr_h(n,:), lags] = xcov(C_err_h-sub_h, V_err_h-sub_h, 30);
    
%     [V1CA1ecorr(n,3), peVC(n,3)] = corr(sign(V_err_h).*(abs(V_err_h)>thres), sign(C_err_h).*(abs(C_err_h))>thres);
%     [delay_ecorr_h(n,:), elags] = xcov(sign(V_err_h).*(abs(V_err_h)>thres), sign(C_err_h).*(abs(C_err_h)>thres), 30);
    [V1CA1ecorr(n,3), peVC(n,3)] = corr((abs(V_err_h)>thres), (abs(C_err_h)>thres));
    [delay_ecorr_h(n,:), elags] = xcov((abs(V_err_h)>thres), (abs(C_err_h)>thres), 30);
%     [V1CA1ecorr(n,3), peVC(n,3)] = corr((abs(V_err_h)>thres), (abs(C_err_h(randperm(length(C_err_h)))))>thres);
%     [delay_ecorr_h(n,:), elags] = xcov((abs(V_err_h)>thres), (abs(C_err_h)>thres), 30);
    
end

% V1CA1ecorr = V1CA1corr;
% peVC = pVC;
% delay_ecorr_h = delay_corr_h;
% delay_ecorr_n = delay_corr_n;
% delay_ecorr_l = delay_corr_l;
% lags = elags;

figure(10);
subplot(231)
temp = pVC(:,1)<1 & pVC(:,3)<1;
plot(V1CA1ecorr(temp,1), V1CA1ecorr(temp,3), 'ko');
axis tight; 
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
xlims = xlim; ylims = ylim;
axis equal; axis square;
axis([min(V1CA1ecorr(:))-0.02 max(V1CA1ecorr(:))+0.02 ...
    min(V1CA1ecorr(:))-0.02 max(V1CA1ecorr(:))+0.02]);
line(xlim, ylim, 'color','k','linestyle', '--')
title('V1-CA1 error ecorrelation');
xlabel('Low Contrast')
ylabel('High Contrast')

subplot(232)
temp = pVC(:,2)<1 & pVC(:,3)<1;
plot(V1CA1ecorr(temp,2), V1CA1ecorr(temp,3), 'ko');
axis tight; 
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
xlims = xlim; ylims = ylim;
axis equal; axis square;
axis([min(V1CA1ecorr(:))-0.02 max(V1CA1ecorr(:))+0.02 ...
    min(V1CA1ecorr(:))-0.02 max(V1CA1ecorr(:))+0.02]);
line(xlim, ylim, 'color','k','linestyle', '--')
title('V1-CA1 error ecorrelation');
xlabel('Medium Contrast')
ylabel('High Contrast')

subplot(234)
temp = pVC(:,1)<1 & pVC(:,3)<1;
plot(V1CA1ecorr(temp,1), V1CA1ecorr(temp,2), 'ko');
axis tight; 
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
xlims = xlim; ylims = ylim;
axis equal; axis square;
axis([min(V1CA1ecorr(:))-0.02 max(V1CA1ecorr(:))+0.02 ...
    min(V1CA1ecorr(:))-0.02 max(V1CA1ecorr(:))+0.02]);
line(xlim, ylim, 'color','k','linestyle', '--')
title('V1-CA1 error ecorrelation');
xlabel('Low Contrast')
ylabel('Medium Contrast')

subplot(233)
hold off;
errorbar([1 2 3], nanmean(V1CA1ecorr), nansem(V1CA1ecorr),'k');
hold on;
bar([1 2 3],nanmean(V1CA1ecorr),'k')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
set(gca, 'YLim', [0 max(V1CA1ecorr(:))]);
title('V1-CA1 error ecorrelation');
set(gca, 'XTick',[1 2 3], 'XTickLabel',[{'Low'}, {'Med'}, {'High'}])

subplot(236)
hold off;
errorarea_as(lags*16.67, nanmean(delay_ecorr_l), nansem(delay_ecorr_l),'b')
hold on;
errorarea_as(lags*16.67, nanmean(delay_ecorr_n), nansem(delay_ecorr_n),'k')
errorarea_as(lags*16.67, nanmean(delay_ecorr_h), nansem(delay_ecorr_h),'r')
axis tight
line([0 0],ylim)
