function [CCG, lags] = spkCCG(spkTimes1, spkTimes2, maxLag)

if nargin<3
    maxLag = 100;
end
if nargin<2
    spkTimes2 = spkTimes1;
    autoC = 1;
elseif length(spkTimes2) == 1
    maxLag = spkTimes2;
    spkTimes2 = spkTimes1;
    autoC = 1;
else
    autoC = 0;
end

% Calculate msec spike times
spkTimes1 = ceil(spkTimes1*1000);
spkTimes2 = ceil(spkTimes2*1000);

% Make spikeTrain arrays
maxTime = max([spkTimes1(:)' spkTimes2(:)']);
spkTrain1 = zeros(1,maxTime);
spkTrain2 = zeros(1,maxTime);
spkTrain1(spkTimes1) = 1;
spkTrain2(spkTimes2) = 1;

% Calculate the cross-correlation
t1 = tic;
[CCG, lags] = xcorr(spkTrain1, spkTrain2, maxLag, 'coeff');
t = toc(t1);
% if autoC
    CCG(abs(lags)<=1) = nan;
% end
disp(['That took ' num2str(t) ' seconds.']);