function [Posterior_all] = plot_some_decoding_byTheta_new5(animal, iseries,...
    exp_list, type, only_correct, smthWin, box_filt, quickProcess, shank_list, area)

%% Parameters that can be changed
test_normal = 1; % Learning/training on normal condition
plot_all = 0;
numSampleBins = 250;
gen_smth_win = 16;

display(['Smoothing window = ' num2str(smthWin)]);
if nargin<4
    type = 'contrast';
end
if nargin<5
    only_correct = false;
end
if nargin<6
    smthWin = 250;
end
if nargin<7
    box_filt = 0;
end
if nargin<8
    quickProcess = 0;
end
if nargin<9
    shank_list = false;
end
if nargin<10
    area = 'CA1';
end

if length(num2str(animal))<=5
    animal = ['M130' num2str(animal) '_BALL'];
end

if strcmp(animal,'M130918_BALL') && (nargin<2 || isempty(exp_list))
    exp_list = 103:105;
elseif strcmp(animal,'M130920_BALL') && (nargin<2 || isempty(exp_list))
    exp_list = 102:103;
end



switch animal
    case 'M130918_BALL'
        iseries = 1030;
        switch area
            case 'V1'
                cell_list = 190:321;
            case 'CA1'
                cell_list = 1:189;
        end
    case 'M130920_BALL'
        iseries = 1025;
        switch area
            case 'V1'
                cell_list = 153:249;
            case 'CA1'
                cell_list = 1:152;
        end
    case 'M130703_BALL'
        iseries = 809;
        cell_list = [];
    otherwise
        if nargin<2
            iseries = inputdlg('Enter the Series:','0');
            iseries = str2num(iseries{1});
        end
        cell_list = [];
end
% delayT = inputdlg('Enter delay frames:','s');
% delayT = str2num(delayT{1});
delayT = 0;

% SetDirs;
% es = VRLoadMultipleExpts(animal, iseries, exp_list,'SPIKES',[],shank_list);
% if ~isempty(cell_list)
%     es.spikeTrain = es.spikeTrain(:,cell_list);
% end
% es.spikeTrain = circshift(es.spikeTrain,[-delayT 0]);
%% Load the spiking and theta data
% SetDirs;
es = [];
% VRLoadMultipleExpts(animal, iseries, exp_list,'SPIKES_THETA',[18 22],shank_list);
% if ~isempty(cell_list)
%     es.spikeTrain = es.spikeTrain(:,cell_list);
% end
if (iseries)==531
    load(['Data' filesep 'es_' animal '_' num2str(iseries) '_thetaBinsA']);
else
    load(['Data' filesep 'es_' animal '_' num2str(iseries) '_thetaBinsB']);
end

switch iseries
    case 1025
        display('Warning!!! Removing interneurons on 1025');
        interneurons = zeros(1,size(es.spikeTrain,2));
        interneurons([17, 18,  9,11,12, 44, 45, 41, 78, 102, 104, 129, 113, 140, 145]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
        
    case 1030
        display('Warning!!! Removing interneurons on 1030');
        interneurons = zeros(1,size(es.spikeTrain,2));
        interneurons([15, 20, 30, 39, 71, 132]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
        
    case 602
        display('Warning!!! Removing interneurons on 602');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([5 8 10 18 24 29 35 38 40 52 54 59 60 61 77 90 92 93 107 109 127]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    case 604
        display('Warning!!! Removing interneurons on 604');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([8 22 25 26 19 31 40 52 53 66 91 97 109]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    case 531
        display('Warning!!! Removing interneurons on 531');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([11 15 16 44 53 86 94 98]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    case 603
        display('Warning!!! Removing interneurons on 603');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([3 8 9 15 57 63 64 66 74 116 131 138 171 177 182 185 188]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    case 601
        display('Warning!!! Removing interneurons on 601');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([4 9 10 11 12 20 24 29 47 50 52 60 64 71 84]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
    case 530
        display('Warning!!! Removing interneurons on 530');
        interneurons = zeros(1,length(es.spikeIDs));
        interneurons([4 8 36 39 56 62 72 78 84 95 98 101 99 102 103 105 131 155 163 167 207 217 222 229 231 266 270]) = 1;
        es.spikeIDs = es.spikeIDs(~interneurons);
        es.spikeTrain = es.spikeTrain(:,~interneurons);
end

% spdDec = bayesDecoder;
% spdDec.CVO = [];
% temp1 = es.smthBallSpd;
% % temp1(temp1<5) = nan;
% [spdDec,a] = spdDec.trainDecoder(temp1, es.spikeTrain,0);
% a = [a' nan*ones(1,length(es.smthBallSpd)-length(a))]';
% es.spikeTrain = a;
%
% es.orig.spikeTrain = es.orig.smthBallSpd;

display(['Processing: ' animal '_' num2str(iseries)]);
%% Getting conditions
Posterior_all.animal = animal;
Posterior_all.iseries = iseries;
Posterior_all.exp_list = exp_list;

switch type
    case 'contrast'
        base = es.traj>1 & es.contrast~=0 & ~isnan(es.traj) & round(100*es.gain)/100==1 & round(100*es.roomLength)/100==1 & es.smthBallSpd>5 & es.trajspeed>=0;
        t      = round(es.contrast*100)/100==0.6 & base;
        t_low  = round(es.contrast*100)/100==0.18 & base;
        t_high = round(es.contrast*100)/100==0.72 & base;
        t_gray = es.traj~=0 & es.contrast==0 & ~isnan(es.traj) & es.smthBallSpd>5;
        if strcmp(animal,'M130703_BALL');
            t      = es.contrast==0.75 & base;
            t_low  = es.contrast<0.75 & base;
            t_high = es.contrast>0.75 & base;
        end
        
        base_orig = es.orig.traj~=0 & es.orig.contrast~=0 & ~isnan(es.orig.traj) & es.orig.gain==1 & es.orig.roomLength==1 & es.orig.smthBallSpd>5 & es.orig.trajspeed>=0;
        t_orig      = es.orig.contrast==0.6 & base_orig;
        t_low_orig  = es.orig.contrast==0.18 & base_orig;
        t_high_orig = es.orig.contrast==0.72 & base_orig;
        t_gray_orig = es.orig.traj~=0 & es.orig.contrast==0 & ~isnan(es.orig.traj) & es.orig.smthBallSpd>5;
        if strcmp(animal,'M130703_BALL');
            t_orig      = es.orig.contrast==0.75 & base;
            t_low_orig  = es.orig.contrast<0.75 & base;
            t_high_orig = es.orig.contrast>0.75 & base;
        end
    case 'gain'
        base = es.traj>1 & round(es.contrast*100)/100==0.6 & ~isnan(es.traj) & round(100*es.gain)/100==1 & es.smthBallSpd>5 & es.trajspeed>=0;
        t      = es.gain==1 & base;
        t_low  = es.gain<1 & base;
        t_high = es.gain>1 & base;
    case 'roomlength'
        base = es.traj>1 & round(es.contrast*100)/100==0.6 & ~isnan(es.traj) & round(100*es.gain)/100==1 & es.smthBallSpd>5 & es.trajspeed>=0;
        t      = es.roomLength==1 & base;
        t_low  = es.roomLength<1 & base;
        t_high = es.roomLength>1 & base;
end

spkRate = zeros(size(es.spikeTrain));
es.spikeTrain(es.traj==0,:) = nan; % to see if there is an effect of the gray screen.

spkRate_orig = zeros(size(es.orig.spikeTrain));
es.orig.spikeTrain(es.orig.traj==0,:) = nan; % to see if there is an effect of the gray screen.
smthWin2 = 125;

for icell = 1:size(es.spikeTrain,2);
    if box_filt
        spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), 7, smthWin,'box');
    else
        spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), 7, smthWin); %*(smthWin*60./1000);
    end
    if box_filt
        spkRate_orig(:,icell) = smthInTime(es.orig.spikeTrain(:,icell), 60, smthWin2,'box');
    else
        spkRate_orig(:,icell) = smthInTime(es.orig.spikeTrain(:,icell), 60, smthWin2); %*(smthWin*60./1000);
    end
end
%% Decoding all the conditions

dec = bayesDecoder;
dec.numBins = numSampleBins; dec.fixedSmth = gen_smth_win;

% if ~test_normal
dec_l = bayesDecoder;%dec;
dec_l.numBins = numSampleBins; dec_l.fixedSmth = gen_smth_win;
dec_h = bayesDecoder;%dec;
dec_h.numBins = numSampleBins; dec_h.fixedSmth = gen_smth_win;

dec_l_d = bayesDecoder;%dec;
dec_l_d.numBins = numSampleBins; dec_l_d.fixedSmth = gen_smth_win/2;
dec_h_d = bayesDecoder;%dec;
dec_h_d.numBins = numSampleBins; dec_h_d.fixedSmth = gen_smth_win/2;
dec_n_d = bayesDecoder;%dec;
dec_n_d.numBins = numSampleBins; dec_n_d.fixedSmth = gen_smth_win/2;

% end
if only_correct
    dec2 = bayesDecoder;%dec;
    dec2.numBins = numSampleBins; dec2.fixedSmth = gen_smth_win;
    [dec, ~, ~, ~] = ...
        ...%          dec.trainDecoder( es.orig.traj(t_orig & es.orig.outcome==2), spkRate_orig(t_orig & es.orig.outcome==2,:), 0);
        dec.trainDecoder( es.traj(t & es.outcome==2), spkRate(t & es.outcome==2,:), 0);
    [~,  Posterior_norm, ~, nonNormPosterior_norm] = dec.predictBayesDecoder(spkRate(t,:), 0,'mean');
    [dec2, ~, X_norm, ~] = ...
        dec2.trainDecoder( es.traj(t), spkRate(t,:), 0);
    %     Posterior_norm = (2.^Posterior_norm)/dec.numBins;
    %     Posterior_norm(Posterior_norm<-5.6) = -5.6;
else
    dec2 = bayesDecoder;%dec;
    %     dec2.numBins = numSampleBins; dec2.fixedSmth = gen_smth_win;
    %     [dec, ~, ~, ~] = ...
    [dec, ~, X_norm, Posterior_norm, ~, nonNormPosterior_norm] = ...
        ...dec.trainDecoder( es.orig.traj(t_orig), spkRate_orig(t_orig,:), 0);
        dec.trainDecoder( es.traj(t), spkRate(t,:), 0);
    %     [~,  Posterior_norm] = dec.predictBayesDecoder(spkRate(t,:), 0,'mean');
    %     [dec2, ~, X_norm, ~] = ...
    %         dec2.trainDecoder( es.traj(t), spkRate(t,:), 0);
    %         [dec, ~, X_norm, Posterior_norm] = ...
    %         dec.trainDecoder( es.traj(t), spkRate(t,:), 0);
    %     Posterior_norm(Posterior_norm<-5.6) = -5.6;
    %     Posterior_norm = (2.^Posterior_norm)/dec.numBins;
end

if ~test_normal
    [dec_l, ~, X_low, Posterior_low, ~, nonNormPosterior_low] = ...
        dec_l.trainDecoder( es.traj(t_low), spkRate(t_low,:), 0);
    [dec_h, ~, X_high, Posterior_high,  ~, nonNormPosterior_high] = ...
        dec_h.trainDecoder( es.traj(t_high), spkRate(t_high,:), 0);
    %     Posterior_low(Posterior_low<-5.6) = -5.6;
    %     Posterior_high(Posterior_high<-5.6) = -5.6;
    Posterior_high_orig =[];
    Posterior_low_orig  =[];
    %     Posterior_low  = (2.^Posterior_low)/dec.numBins;
    %     Posterior_high = (2.^Posterior_high)/dec.numBins;
else
    [dec_l, ML_low_orig, X_low_orig, Posterior_low_orig, ~, nonNormPosterior_low_orig] = ...
        dec_l.trainDecoder( es.traj(t_low), spkRate(t_low,:), 0);
    [dec_h, ML_high_orig, X_high_orig, Posterior_high_orig, ~, nonNormPosterior_high_orig] = ...
        dec_h.trainDecoder( es.traj(t_high), spkRate(t_high,:), 0);
    
    %     X_high_orig = X_high_orig(1:length(ML_high_orig));
    %     X_low_orig  = X_low_orig(1:length(ML_low_orig));
    
    [ML_low,  Posterior_low, ~, nonNormPosterior_low] = dec.predictBayesDecoder(spkRate(t_low,:), 0,'mean');
    [~,  Posterior_gray, ~, nonNormPosterior_gray] = dec.predictBayesDecoder(spkRate(t_gray,:), 0,'mean');
    [ML_high,  Posterior_high, ~, nonNormPosterior_high] = dec.predictBayesDecoder(spkRate(t_high,:), 0,'mean');
    [ML_norm] = dec.predictBayesDecoder(spkRate(t,:), 0,'mean');
    %     Posterior_low(Posterior_low<-5.6) = -5.6;
    %     Posterior_high(Posterior_high<-5.6) = -5.6;
    %     Posterior_gray(Posterior_gray<-5.6) = -5.6;
    %     Posterior_low  = (2.^Posterior_low)/dec.numBins;
    %     Posterior_high = (2.^Posterior_high)/dec.numBins;
    %     Posterior_gray = (2.^Posterior_gray)/dec.numBins;
    %% Train the decoder on decoded position
    %     dec_l_d = dec_l; % decoder (low contrast trained on decoding position)
    %     dec_h_d = dec_h; % decoder (high contrast trained on decoding position)
    %     dec_n_d = dec;
    %     quickProcess = 1;
    if quickProcess
        dec_l_d = bayesDecoder;
        dec_l_d.numBins = numSampleBins; dec_l_d.fixedSmth = gen_smth_win/2;
        dec_h_d = bayesDecoder;
        dec_h_d.numBins = numSampleBins; dec_h_d.fixedSmth = gen_smth_win/2;
        dec_n_d = bayesDecoder;
        dec_h_d.numBins = numSampleBins; dec_h_d.fixedSmth = gen_smth_win/2;
        
        [dec_l_d, ~, X_low_orig_d, Posterior_low_orig_d] = ...
            dec_l_d.trainDecoder( ML_low', spkRate(t_low,:), 0);
        [dec_h_d, ~, X_high_orig_d, Posterior_high_orig_d] = ...
            dec_h_d.trainDecoder( ML_high', spkRate(t_high,:), 0);
        [dec_n_d, ~, X_norm_orig_d, Posterior_norm_orig_d] = ...
            dec_n_d.trainDecoder( ML_norm', spkRate(t,:), 0);
        
        %         Posterior_low_orig_d(Posterior_low_orig_d<-5.6) = -5.6;
        %         Posterior_high_orig_d(Posterior_high_orig_d<-5.6) = -5.6;
        %         Posterior_norm_orig_d(Posterior_norm_orig_d<-5.6) = -5.6;
        Posterior_low_orig_d  = (2.^Posterior_low_orig_d)/dec.numBins;
        Posterior_high_orig_d = (2.^Posterior_high_orig_d)/dec.numBins;
        Posterior_norm_orig_d = (2.^Posterior_norm_orig_d)/dec.numBins;
    else
        %         tmpDec_low = dec_l;
        %         tmpDec_high = dec_h;
        
        tmpMap = oneDimMap;
        tmpMap.kfold = dec.kfold;
        tmpMap.bins  = dec.bins;
        tmpMap.numBins = dec.numBins;
        
        %         dec_l_d.bestModel_orig = dec_l_d.model.bestModel;
        %         dec_h_d.bestModel_orig = dec_l_d.model.bestModel;
        %         dec_n_d.bestModel_orig = dec_l_d.model.bestModel;
        
        clear temp_dec_*
        for icell = 1:size(spkRate,2)
            display(['Processing cell: ' num2str(icell) '/' num2str(size(spkRate,2)) ' series:' num2str(iseries)]);
            tmpSpkRate = spkRate;
            tmpSpkRate(:,icell) = 0;
            temp_dec = dec;
            temp_dec_l = dec_l;
            temp_dec_h = dec_h;
            %             temp_dec.model.meanModel(icell,:) = 0;
            %             temp_dec_l.model.meanModel(icell,:) = 0;
            %             temp_dec_h.model.meanModel(icell,:) = 0;
            %
            %             temp_dec.model.bestModel(icell,:) = 0;
            %             temp_dec_l.model.bestModel(icell,:) = 0;
            %             temp_dec_h.model.bestModel(icell,:) = 0;
            
            [tmpML_norm, P] = temp_dec.predictBayesDecoder(tmpSpkRate(t,:), 0,'mean');
            %             [tmpML_low ] = dec.predictBayesDecoder(tmpSpkRate(t_low,:), 0,'best');
            %             [tmpML_high] = dec.predictBayesDecoder(tmpSpkRate(t_high,:), 0,'best');
            [tmpML_low_norm ] = temp_dec.predictBayesDecoder(tmpSpkRate(t_low,:), 0,'mean');
            [tmpML_high_norm] = temp_dec.predictBayesDecoder(tmpSpkRate(t_high,:), 0,'mean');
            
            [tmpML_low_orig ] = temp_dec_l.predictBayesDecoder(tmpSpkRate(t_low,:), 0,'mean');
            [tmpML_high_orig] = temp_dec_h.predictBayesDecoder(tmpSpkRate(t_high,:), 0,'mean');
            
            % Choosing Best decoder
            if mean((ML_low_orig - X_low_orig(1:length(ML_low_orig))).^2) < mean((ML_low - X_low_orig').^2)
                tmpML_low = tmpML_low_orig;
            else
                tmpML_low = tmpML_low_norm;
            end
            if mean((ML_high_orig - X_high_orig(1:length(ML_high_orig))).^2) < mean((ML_high - X_high_orig').^2)
                tmpML_high = tmpML_high_orig;
            else
                tmpML_high = tmpML_high_norm;
            end
            
            while length(tmpML_low)<sum(t_low)
                tmpML_low = [tmpML_low nan];
            end
            while length(tmpML_high)<sum(t_high)
                tmpML_high = [tmpML_high nan];
            end
            while length(tmpML_norm)<sum(t)
                tmpML_norm = [tmpML_norm nan];
            end
            
            while length(tmpML_low)<sum(t_low)
                tmpML_low = [tmpML_low nan];
            end
            while length(tmpML_high)<sum(t_high)
                tmpML_high = [tmpML_high nan];
            end
            while length(tmpML_norm)<sum(t)
                tmpML_norm = [tmpML_norm nan];
            end
            Posterior_all.MAP.dec.low(icell,:)  = tmpML_low_norm;
            Posterior_all.MAP.dec.norm(icell,:) = tmpML_norm;
            Posterior_all.MAP.dec.high(icell,:) = tmpML_high_norm;
            
            tmpMap_low = oneDimMap;
            tmpMap_low.kfold = dec.kfold;
            tmpMap_low.bins  = dec.bins;
            tmpMap_low.numBins = dec.numBins;
            
            tmpMap_high = oneDimMap;
            tmpMap_high.kfold = dec.kfold;
            tmpMap_high.bins  = dec.bins;
            tmpMap_high.numBins = dec.numBins;
            
            tmpMap_norm = oneDimMap;
            tmpMap_norm.kfold = dec.kfold;
            tmpMap_norm.bins  = dec.bins;
            tmpMap_norm.numBins = dec.numBins;
            
            tmpMap_low_orig = oneDimMap;
            tmpMap_low_orig.kfold = dec.kfold;
            tmpMap_low_orig.bins  = dec.bins;
            tmpMap_low_orig.numBins = dec.numBins;
            
            tmpMap_high_orig = oneDimMap;
            tmpMap_high_orig.kfold = dec.kfold;
            tmpMap_high_orig.bins  = dec.bins;
            tmpMap_high_orig.numBins = dec.numBins;
            
            tmpMap_norm_orig = oneDimMap;
            tmpMap_norm_orig.kfold = dec.kfold;
            tmpMap_norm_orig.bins  = dec.bins;
            tmpMap_norm_orig.numBins = dec.numBins;
            
            %This limits the decoder to half the track
            rewlim = 50*numSampleBins/100;
            
            % Calculating decoded maps
            %             tmpMap_low = tmpMap;
            %             tmpMap_high = tmpMap;
            tmpMap_low.CVO = [];
            temp = zeros(size(t_low));
            temp(t_low) = X_low_orig;
            temp2 = zeros(size(t_low));
            temp2(t_low) = tmpML_low;
            [tmpMap_low] = tmpMap_low.trainSpikeMap(tmpML_low(X_low_orig'<rewlim & tmpML_low<rewlim)', spkRate(t_low & temp<rewlim & temp2<rewlim,icell), 0, gen_smth_win/2);
            
            tmpMap_high.CVO = [];
            temp = zeros(size(t_high));
            temp(t_high) = X_high_orig;
            %             [tmpMap_high] = tmpMap_high.trainSpikeMap(tmpML_high(X_high_orig<rewlim)', spkRate(t_high & temp<rewlim,icell), 0);
            temp2 = zeros(size(t_high));
            temp2(t_high) = tmpML_high;
            [tmpMap_high] = tmpMap_high.trainSpikeMap(tmpML_high(X_high_orig'<rewlim & tmpML_high<rewlim)', spkRate(t_high & temp<rewlim & temp2<rewlim,icell), 0, gen_smth_win/2);
            
            tmpMap_norm.CVO = [];
            temp = zeros(size(t));
            temp(t) = X_norm;
            %             [tmpMap_norm] = tmpMap_norm.trainSpikeMap(tmpML_norm(X_norm<rewlim)', spkRate(t& temp<rewlim,icell), 0);
            temp2 = zeros(size(t));
            temp2(t) = tmpML_norm;
            [tmpMap_norm] = tmpMap_norm.trainSpikeMap(tmpML_norm(X_norm'<rewlim & tmpML_norm<rewlim)', spkRate(t & temp<rewlim & temp2<rewlim,icell), 0, gen_smth_win/2);
            
            % Recalculating original maps
            tmpMap_low_orig.CVO = [];
            temp = zeros(size(t_low));
            temp(t_low) = X_low_orig;
            [tmpMap_low_orig] = tmpMap_low_orig.trainSpikeMap(X_low_orig(X_low_orig<rewlim)', spkRate(t_low & temp<rewlim,icell), 0, gen_smth_win/2);
            
            tmpMap_high_orig.CVO = [];
            temp = zeros(size(t_high));
            temp(t_high) = X_high_orig;
            [tmpMap_high_orig] = tmpMap_high_orig.trainSpikeMap(X_high_orig(X_high_orig<rewlim)', spkRate(t_high & temp<rewlim,icell), 0, gen_smth_win/2);
            
            tmpMap_norm_orig.CVO = [];
            temp = zeros(size(t));
            temp(t) = X_norm;
            [tmpMap_norm_orig] = tmpMap_norm_orig.trainSpikeMap(X_norm(X_norm<rewlim)', spkRate(t& temp<rewlim,icell), 0, gen_smth_win/2);
            
            dec_l_d.model.EV_orig(:,icell) = tmpMap_low_orig.model.EV; %EV_low; %
            dec_h_d.model.EV_orig(:,icell) = tmpMap_high_orig.model.EV;%EV_high; %
            dec_n_d.model.EV_orig(:,icell) = tmpMap_norm_orig.model.EV;%EV_norm; %
            
            [~, low_idx_orig] = max(tmpMap_low_orig.model.EV);
            [~, high_idx_orig] = max(tmpMap_high_orig.model.EV);
            [~, norm_idx_orig] = max(tmpMap_norm_orig.model.EV);
            
            dec_l_d.model.EV(:,icell) = tmpMap_low.model.EV; %EV_low; %
            dec_h_d.model.EV(:,icell) = tmpMap_high.model.EV;%EV_high; %
            dec_n_d.model.EV(:,icell) = tmpMap_norm.model.EV;%EV_norm; %
            
            [~, low_idx] = max(tmpMap_low.model.EV);
            [~, high_idx] = max(tmpMap_high.model.EV);
            [~, norm_idx] = max(tmpMap_norm.model.EV);
            
            dec_l_d.model.bestModel(icell,:) = tmpMap_low.model.tuning.respModel(low_idx,:);
            dec_h_d.model.bestModel(icell,:) = tmpMap_high.model.tuning.respModel(high_idx,:);
            dec_n_d.model.bestModel(icell,:) = tmpMap_norm.model.tuning.respModel(norm_idx,:);
            
            dec_l_d.model.bestModel_orig(icell,:) = tmpMap_low_orig.model.tuning.respModel(low_idx,:);
            dec_h_d.model.bestModel_orig(icell,:) = tmpMap_high_orig.model.tuning.respModel(high_idx,:);
            dec_n_d.model.bestModel_orig(icell,:) = tmpMap_norm_orig.model.tuning.respModel(norm_idx,:);
            
            for tmpIdx = 1:tmpMap.kfold
                dec_l_d.model.trained(tmpIdx).respModel_orig(icell,:) = tmpMap_low_orig.model.tuning.respModel(tmpIdx,:);
                dec_l_d.model.trained(tmpIdx).respModel(icell,:) = tmpMap_low.model.tuning.respModel(tmpIdx,:);
                tempModel_low(tmpIdx,icell,:) = dec_l_d.model.trained(tmpIdx).respModel(icell,:);
                tempModel_low_orig(tmpIdx,icell,:) = dec_l_d.model.trained(tmpIdx).respModel_orig(icell,:);
                
                dec_h_d.model.trained(tmpIdx).respModel_orig(icell,:) = tmpMap_high_orig.model.tuning.respModel(tmpIdx,:);
                dec_h_d.model.trained(tmpIdx).respModel(icell,:) = tmpMap_high.model.tuning.respModel(tmpIdx,:);
                tempModel_high(tmpIdx,icell,:) = dec_h_d.model.trained(tmpIdx).respModel(icell,:);
                tempModel_high_orig(tmpIdx,icell,:) = dec_h_d.model.trained(tmpIdx).respModel_orig(icell,:);
                
                dec_n_d.model.trained(tmpIdx).respModel_orig(icell,:) = tmpMap_norm_orig.model.tuning.respModel(tmpIdx,:);
                dec_n_d.model.trained(tmpIdx).respModel(icell,:) = tmpMap_norm.model.tuning.respModel(tmpIdx,:);
                tempModel_norm(tmpIdx,icell,:) = dec_n_d.model.trained(tmpIdx).respModel(icell,:);
                tempModel_norm_orig(tmpIdx,icell,:) = dec_n_d.model.trained(tmpIdx).respModel_orig(icell,:);
            end
            clear tmpSpkRate tmpMap_low tmpMap_high tmpMap_norm;
            drawnow
        end
        dec_l_d.model.meanModel = squeeze(nanmean(tempModel_low,1));
        dec_h_d.model.meanModel = squeeze(nanmean(tempModel_high,1));
        dec_n_d.model.meanModel = squeeze(nanmean(tempModel_norm,1));
        
        dec_l_d.model.meanModel_orig = squeeze(nanmean(tempModel_low_orig,1));
        dec_h_d.model.meanModel_orig = squeeze(nanmean(tempModel_high_orig,1));
        dec_n_d.model.meanModel_orig = squeeze(nanmean(tempModel_norm_orig,1));
    end
end
%%
if test_normal
    [X_low, bins_low] = normalise1var(es.traj(t_low), dec.numBins);
    [X_high,bins_high] = normalise1var(es.traj(t_high), dec.numBins);
    X_high_orig = X_high(1:size(Posterior_high_orig,1));
    X_low_orig  = X_low(1:size(Posterior_low_orig,1));
else
    X_high = X_high(1:size(Posterior_high,1));
    bins_low = dec_l.bins;
    X_low  = X_low(1:size(Posterior_low,1));
    bins_high = dec_h.bins;
end
X_norm = X_norm(1:size(Posterior_norm,1));

timePs = find(t);
timePs = timePs(1:size(Posterior_norm,1));
t = false(size(t));
t(timePs) = true;

timePs = find(t_low);
timePs = timePs(1:size(Posterior_low,1));
t_low = false(size(t_low));
t_low(timePs) = true;

timePs = find(t_high);
timePs = timePs(1:size(Posterior_high,1));
t_high = false(size(t_high));
t_high(timePs) = true;
%% Defining the outputs
Posterior_all.Posterior_norm = Posterior_norm;
Posterior_all.Posterior_high = Posterior_high;
Posterior_all.Posterior_low  = Posterior_low ;

Posterior_all.nonNormPosterior_norm = nonNormPosterior_norm;
Posterior_all.nonNormPosterior_high = nonNormPosterior_high;
Posterior_all.nonNormPosterior_low  = nonNormPosterior_low ;
Posterior_all.nonNormPosterior_gray = nonNormPosterior_gray;


Posterior_all.Posterior_gray  = Posterior_gray;
Posterior_all.only_correct = only_correct;

Posterior_all.X_norm = X_norm;
Posterior_all.X_high = X_high;
Posterior_all.X_low  = X_low;

Posterior_all.data      = es;
Posterior_all.decoder   = dec;
Posterior_all.decoder_low   = dec_l;
Posterior_all.decoder_high   = dec_h;

Posterior_all.decoder_low_decPos   = dec_l_d;
Posterior_all.decoder_high_decPos   = dec_h_d;
Posterior_all.decoder_norm_decPos   = dec_n_d;

Posterior_all.t_norm    = t;
Posterior_all.t_low     = t_low;
Posterior_all.t_high    = t_high;

Posterior_all.meanrates.low = mean(es.spikeTrain(t_low & es.traj<60,:),1);
Posterior_all.meanrates.norm = mean(es.spikeTrain(t & es.traj<60,:),1);
Posterior_all.meanrates.high = mean(es.spikeTrain(t_high & es.traj<60,:),1);

Posterior_all.Posterior_high_orig = Posterior_high_orig;
Posterior_all.Posterior_low_orig  = Posterior_low_orig;
Posterior_all.X_high_orig = X_high_orig;
Posterior_all.X_low_orig  = X_low_orig;
%% Finding the trials as correct, early, miss and late
trialEnds = [find(diff(es.trialID)>=1)-15];
temp = zeros(size(es.traj));
temp(trialEnds) = 1;
trialEnds = temp;

outcome.complete        = (es.trajPercent(trialEnds>0)>80) | es.outcome(trialEnds>0)~=1;
outcome.correctTrials   = es.trialID(trialEnds & es.outcome==2);
outcome.earlyTrials     = es.trialID(trialEnds & (es.outcome==0 & es.trajPercent<80));
outcome.misslateTrials  = es.trialID(trialEnds  & es.outcome==0 & es.trajPercent>80);
outcome.lowContrast     = es.trialID(trialEnds & t_low);
outcome.normContrast    = es.trialID(trialEnds & t);
outcome.highContrast    = es.trialID(trialEnds & t_high);

es.outcome(es.outcome==1) = NaN;

% Splitting miss and late trials
outcome.missTrials = [];
outcome.lateTrials = [];
for tIdx = outcome.misslateTrials'
    % Does the animal lick after the reward position?
    if sum(es.lick((es.trialID==tIdx) & (es.trajPercent>70)))>0
        outcome.lateTrials = [outcome.lateTrials tIdx];
    else
        outcome.missTrials = [outcome.missTrials tIdx];
    end
end

% Early
for tIdx = [outcome.earlyTrials']
    es.outcome(es.trialID==tIdx) = 1;
end

% Correct
for tIdx = [outcome.correctTrials']
    es.outcome(es.trialID==tIdx) = 2;
end
% Late
for tIdx = [outcome.lateTrials]
    es.outcome(es.trialID==tIdx) = 3;
end
% Miss
for tIdx = [outcome.missTrials]
    es.outcome(es.trialID==tIdx) = 4;
end
%% Getting the confusion matrixes
for n = 1:dec.numBins
    %
    timePts = (X_norm == n);
    meanPost_norm(:,n) = nanmean(Posterior_norm(timePts, :),1);
    %
    timePts = (X_high == n);
    meanPost_high(:,n) = nanmean(Posterior_high(timePts, :),1);
    %
    timePts = (X_low == n);
    meanPost_low(:,n) = nanmean(Posterior_low(timePts, :),1);
    %
    timePts = (X_high_orig == n);
    meanPost_high_orig(:,n) = nanmean(Posterior_high_orig(timePts, :),1);
    %
    timePts = (X_low_orig == n);
    meanPost_low_orig(:,n) = nanmean(Posterior_low_orig(timePts, :),1);
    %
    timePts = (X_norm == n & (es.outcome(t) ~= 2));
    meanPost_norm_pass(:,n) = nanmean(Posterior_norm(timePts, :),1);
    %
    timePts = (X_norm == n & (es.outcome(t) == 2));
    meanPost_norm_actv(:,n) = nanmean(Posterior_norm(timePts, :),1);
    %
    timePts = (X_low == n & (es.outcome(t_low) == 2));
    meanPost_low_actv(:,n) = nanmean(Posterior_low(timePts, :),1);
    %
    timePts = (X_low == n & (es.outcome(t_low) ~= 2));
    meanPost_low_pass(:,n) = nanmean(Posterior_low(timePts, :),1);
    %
    timePts = (X_high == n & (es.outcome(t_high) ~= 2));
    meanPost_high_pass(:,n) = nanmean(Posterior_high(timePts, :),1);
    %
    timePts = (X_high == n & (es.outcome(t_high) == 2));
    meanPost_high_actv(:,n) = nanmean(Posterior_high(timePts, :),1);
end

outcome.allContrast = [outcome.lowContrast' outcome.normContrast' outcome.highContrast']';
%% Getting the Early, Correct, Late and Miss trials
for n = 1:dec.numBins
    % Early
    % % Low
    timePts1 = (X_low == n & (es.outcome(t_low) == 1));
    meanPost.low.early(:,n) = nanmean(Posterior_low(timePts1, :),1);
    % % Norm
    timePts2 = (X_norm == n & (es.outcome(t) == 1));
    meanPost.norm.early(:,n) = nanmean(Posterior_norm(timePts2, :),1);
    % % High
    timePts3 = (X_high == n & (es.outcome(t_high) == 1));
    meanPost.high.early(:,n) = nanmean(Posterior_high(timePts3, :),1);
    % % All
    meanPost.all.early(:,n) = nanmean([Posterior_low(timePts1, :)' Posterior_norm(timePts2, :)' Posterior_high(timePts3, :)']',1);
    
    % Correct
    % % Low
    timePts1 = (X_low == n & (es.outcome(t_low) == 2));
    meanPost.low.correct(:,n) = nanmean(Posterior_low(timePts1, :),1);
    % % Norm
    timePts2 = (X_norm == n & (es.outcome(t) == 2));
    meanPost.norm.correct(:,n) = nanmean(Posterior_norm(timePts2, :),1);
    % % High
    timePts3 = (X_high == n & (es.outcome(t_high) == 2));
    meanPost.high.correct(:,n) = nanmean(Posterior_high(timePts3, :),1);
    % % All
    meanPost.all.correct(:,n) = nanmean([Posterior_low(timePts1, :)' Posterior_norm(timePts2, :)' Posterior_high(timePts3, :)']',1);
    
    % Late
    % % Low
    timePts1 = (X_low == n & (es.outcome(t_low) == 3));
    meanPost.low.late(:,n) = nanmean(Posterior_low(timePts1, :),1);
    % % Norm
    timePts2 = (X_norm == n & (es.outcome(t) == 3));
    meanPost.norm.late(:,n) = nanmean(Posterior_norm(timePts2, :),1);
    % % High
    timePts3 = (X_high == n & (es.outcome(t_high) == 3));
    meanPost.high.late(:,n) = nanmean(Posterior_high(timePts3, :),1);
    % % All
    meanPost.all.late(:,n) = nanmean([Posterior_low(timePts1, :)' Posterior_norm(timePts2, :)' Posterior_high(timePts3, :)']',1);
    
    % Miss
    % % Low
    timePts1 = (X_low == n & (es.outcome(t_low) == 4));
    meanPost.low.miss(:,n) = nanmean(Posterior_low(timePts1, :),1);
    % % Norm
    timePts2 = (X_norm == n & (es.outcome(t) == 4));
    meanPost.norm.miss(:,n) = nanmean(Posterior_norm(timePts2, :),1);
    % % High
    timePts3 = (X_high == n & (es.outcome(t_high) == 4));
    meanPost.high.miss(:,n) = nanmean(Posterior_high(timePts3, :),1);
    % % All
    meanPost.all.miss(:,n) = nanmean([Posterior_low(timePts1, :)' Posterior_norm(timePts2, :)' Posterior_high(timePts3, :)']',1);
    
end
%%
Posterior_all.outcome = outcome;

Posterior_all.meanPost.norm = meanPost_norm;
Posterior_all.meanPost.high = meanPost_high;
Posterior_all.meanPost.low =  meanPost_low;

Posterior_all.meanPost.high_orig = meanPost_high_orig;
Posterior_all.meanPost.low_orig =  meanPost_low_orig;

Posterior_all.meanPost.norm_pass = meanPost_norm_pass;
Posterior_all.meanPost.high_pass = meanPost_high_pass;
Posterior_all.meanPost.low_pass =  meanPost_low_pass;

Posterior_all.meanPost.norm_actv = meanPost_norm_actv;
Posterior_all.meanPost.high_actv = meanPost_high_actv;
Posterior_all.meanPost.low_actv =  meanPost_low_actv;

Posterior_all.meanPost_new.low.early    =  meanPost.low.early;
Posterior_all.meanPost_new.low.correct  =  meanPost.low.correct;
Posterior_all.meanPost_new.low.late     =  meanPost.low.late;
Posterior_all.meanPost_new.low.miss     =  meanPost.low.miss;

Posterior_all.meanPost_new.norm.early    =  meanPost.norm.early;
Posterior_all.meanPost_new.norm.correct  =  meanPost.norm.correct;
Posterior_all.meanPost_new.norm.late     =  meanPost.norm.late;
Posterior_all.meanPost_new.norm.miss     =  meanPost.norm.miss;

Posterior_all.meanPost_new.high.early    =  meanPost.high.early;
Posterior_all.meanPost_new.high.correct  =  meanPost.high.correct;
Posterior_all.meanPost_new.high.late     =  meanPost.high.late;
Posterior_all.meanPost_new.high.miss     =  meanPost.high.miss;

Posterior_all.meanPost_new.all.early    =  meanPost.all.early;
Posterior_all.meanPost_new.all.correct  =  meanPost.all.correct;
Posterior_all.meanPost_new.all.late     =  meanPost.all.late;
Posterior_all.meanPost_new.all.miss     =  meanPost.all.miss;

bestModel = dec.model.bestModel;
bestModel(find(sum(60*bestModel')<1),:) = [];
[~,maxPos] = max(bestModel');
[~,sort_order] = sort(maxPos);
for n = 1:size(bestModel,1)
    bestModel(n,:) = bestModel(n,:) - min(bestModel(n,:));
    bestModel(n,:) = bestModel(n,:) ./max(bestModel(n,:));
end

Posterior_all.data      = es;

if plot_all
    %% Plotting from here
    figure(5)
    subplot(311)
    hold off;
    imagesc(Posterior_low')
    RedWhiteBlue;
    hold on
    axis xy
    plot(1:sum(t_low),X_low,'k.', 'linewidth',2)
    plot(find(es.lick(t_low)==1),X_low(es.lick(t_low)==1),'mo')
    plot(find(es.lick(t_low)==1),X_low(es.lick(t_low)==1),'g.')
    plot(find(es.lick(t_low)==1),X_low(es.lick(t_low)==1),'mo')
    subplot(312)
    hold off;
    imagesc(Posterior_norm')
    axis xy
    hold on;
    plot(1:sum(t),X_norm,'k.', 'linewidth',2)
    plot(find(es.lick(t)==1),X_norm(es.lick(t)==1),'g.')
    plot(find(es.lick(t)==1),X_norm(es.lick(t)==1),'mo')
    subplot(313)
    hold off
    imagesc(Posterior_high')
    axis xy
    hold on;
    plot(1:sum(t_high),X_high,'k.', 'linewidth',2)
    plot(find(es.lick(t_high)==1),X_high(es.lick(t_high)==1),'g.')
    plot(find(es.lick(t_high)==1),X_high(es.lick(t_high)==1),'mo')
    for n = 1:3
        subplot(3,1,n)
        % set(gca,'CLim',[-5 5])
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    end
    % set(gca,'CLim',[-5 5])
    hcb = colorbar('YTick',[-5 5],'YTickLabel',{'2^5 x', 'Chance', '2^5 x'});
    set(hcb,'YTickMode','manual')
    % colormap(gray);
    %%
    figure(1)
    imagesc(bestModel(sort_order,:)); axis xy;
    
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    ylabel('Cell #')
    xlabel('Position in room')
    %%
    figure(2)
    subplot(311)
    imagesc(bins_low, dec.bins, meanPost_low); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['Low ' type], 'fontsize', 14);
    subplot(312)
    imagesc(dec.bins, dec.bins, meanPost_norm); axis xy; colorbar; axis tight; axis equal; axis tight
    title('Normal', 'fontsize', 14);
    subplot(313)
    imagesc(bins_high, dec.bins, meanPost_high); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['High ' type], 'fontsize', 14);
    RedWhiteBlue;
    for n = 1:3
        subplot(3,1,n)
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
        line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
    end
    for n = 1:3
        subplot(3,1,n)
        xlabel('Original Position')
        ylabel('Decoded Posterior')
        % set(gca,'CLim',[-0.5 0.5])
        hcb = colorbar('YTick',[-1 0 1],'YTickLabel',{'2^-0.5 x', 'Chance', '2^0.5 x'});
        set(hcb,'YTickMode','manual')
    end
    
    %%
    figure(8)
    subplot(321)
    imagesc(dec.bins, dec.bins, meanPost_low_pass); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['Low ' type ': MISS'], 'fontsize', 14);
    RedWhiteBlue;
    subplot(322)
    imagesc(dec.bins, dec.bins, meanPost_low_actv); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['Low ' type ': HIT'], 'fontsize', 14);
    subplot(323)
    imagesc(dec.bins, dec.bins, meanPost_norm_pass); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['Normal: MISS'], 'fontsize', 14);
    subplot(324)
    imagesc(dec.bins, dec.bins, meanPost_norm_actv); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['Normal: HIT'], 'fontsize', 14);
    subplot(325)
    imagesc(dec.bins, dec.bins, meanPost_high_pass); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['High ' type ': MISS'], 'fontsize', 14);
    subplot(326)
    imagesc(dec.bins, dec.bins, meanPost_high_actv); axis xy; colorbar; axis tight; axis equal; axis tight
    title(['High ' type ': HIT'], 'fontsize', 14);
    for n = 1:6
        subplot(3,2,n)
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
        line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
    end
    for n = 1:6
        subplot(3,2,n)
        xlabel('Original Position')
        ylabel('Decoded Posterior')
        % set(gca,'CLim',[-0.5 0.5])
        hcb = colorbar('YTick',[-1 0 1],'YTickLabel',{'2^-0.5 x', 'Chance', '2^0.5 x'});
        set(hcb,'YTickMode','manual')
    end
end
%% Calculating error, accuracy and confidence
t_n = zeros(size(Posterior_all.Posterior_norm));
t_l = zeros(size(Posterior_all.Posterior_low));
t_h = zeros(size(Posterior_all.Posterior_high));
for n = 1:size(t_l,1)
    t_l(n,Posterior_all.X_low(n)) = 1;
end
for n = 1:size(t_h,1)
    t_h(n,Posterior_all.X_high(n)) = 1;
end
%     X_norm = X_norm(1:size(Posterior_all.Posterior_norm,1));
for n = 1:size(t_n,1)
    t_n(n,Posterior_all.X_norm(n)) = 1;
end
Posterior_all.confidence.low = (max(Posterior_all.Posterior_low') - min(Posterior_all.Posterior_low'));
Posterior_all.confidence.norm = (max(Posterior_all.Posterior_norm')- min(Posterior_all.Posterior_norm'));
Posterior_all.confidence.high = (max(Posterior_all.Posterior_high') - min(Posterior_all.Posterior_high'));
Posterior_all.confidence.gray = (max(Posterior_all.Posterior_gray') - min(Posterior_all.Posterior_gray'));

Posterior_all.accuracy.low =  Posterior_all.Posterior_low(t_l>0)./max(Posterior_all.Posterior_low(t_l>0)');
Posterior_all.accuracy.norm = Posterior_all.Posterior_norm(t_n>0)./max(Posterior_all.Posterior_norm(t_n>0)');
Posterior_all.accuracy.high = Posterior_all.Posterior_high(t_h>0)./max(Posterior_all.Posterior_high(t_h>0)');

[~,X_ML_low] = max(Posterior_all.Posterior_low');
[~,X_ML_norm] = max(Posterior_all.Posterior_norm');
[~,X_ML_high] = max(Posterior_all.Posterior_high');
[~,X_ML_gray] = max(Posterior_all.Posterior_gray');

Posterior_all.MAP.low = X_ML_low;
Posterior_all.MAP.norm = X_ML_norm;
Posterior_all.MAP.high = X_ML_high;
Posterior_all.MAP.gray = X_ML_gray;

Posterior_all.error.low  = (abs(X_ML_low - X_low'));
Posterior_all.error.norm = (abs(X_ML_norm - X_norm'));
Posterior_all.error.high = (abs(X_ML_high - X_high'));
%     Posterior_all.accuracy.low = mean(Posterior_all.Posterior_low(t_l>0));
%     Posterior_all.accuracy.norm = mean(Posterior_all.Posterior_norm(t>0));
%     Posterior_all.accuracy.high = mean(Posterior_all.Posterior_high(t_h>0));

%     Posterior_all.width.low = mean(std(Posterior_all.Posterior_low'));
%     Posterior_all.width.norm = mean(std(Posterior_all.Posterior_norm'));
%     Posterior_all.width.high = mean(std(Posterior_all.Posterior_high'));

t_l = zeros(size(Posterior_all.Posterior_low_orig));
t_h = zeros(size(Posterior_all.Posterior_high_orig));
for n = 1:size(t_l,1)
    t_l(n,Posterior_all.X_low_orig(n)) = 1;
end
for n = 1:size(t_h,1)
    t_h(n,Posterior_all.X_high_orig(n)) = 1;
end
Posterior_all.confidence.low_orig = (max(Posterior_all.Posterior_low_orig') - min(Posterior_all.Posterior_low_orig'));
Posterior_all.confidence.high_orig = (max(Posterior_all.Posterior_high_orig') - min(Posterior_all.Posterior_high_orig'));

Posterior_all.accuracy.low_orig = (Posterior_all.Posterior_low_orig(t_l>0))./max(Posterior_all.Posterior_low_orig(t_l>0)');
Posterior_all.accuracy.high_orig = (Posterior_all.Posterior_high_orig(t_h>0))./max(Posterior_all.Posterior_high_orig(t_h>0)');

[~,X_ML_low_orig] = max(Posterior_all.Posterior_low_orig');
[~,X_ML_high_orig] = max(Posterior_all.Posterior_high_orig');

Posterior_all.error.low_orig = (abs(X_ML_low_orig - X_low_orig(1:length(X_ML_low_orig))'));;
Posterior_all.error.high_orig = (abs(X_ML_high_orig - X_high_orig(1:length(X_ML_high_orig))'));;

Posterior_all.width.low_orig = mean(std(Posterior_all.Posterior_low_orig'));
Posterior_all.width.high_orig = mean(std(Posterior_all.Posterior_high_orig'));
%% getting the 45 deg marginals
tmp = Posterior_all.meanPost.low(:,1:32);
[Posterior_all.marginals.low, Posterior_all.marginals.lowX] = get45Marginal(tmp);
tmp = Posterior_all.meanPost.high(:,1:32);
[Posterior_all.marginals.high, Posterior_all.marginals.highX] = get45Marginal(tmp);
tmp = Posterior_all.meanPost.norm(:,1:32);
[Posterior_all.marginals.norm, Posterior_all.marginals.normX] = get45Marginal(tmp);

tmp = Posterior_all.meanPost.low;
[Posterior_all.marginals.alllow, Posterior_all.marginals.alllowX] = get45Marginal(tmp);
tmp = Posterior_all.meanPost.high;
[Posterior_all.marginals.allhigh, Posterior_all.marginals.allhighX] = get45Marginal(tmp);
tmp = Posterior_all.meanPost.norm;
[Posterior_all.marginals.allnorm, Posterior_all.marginals.allnormX] = get45Marginal(tmp);
%%
if plot_all
    figure;
    plot(Posterior_all.marginals.lowX*2*sqrt(2),Posterior_all.marginals.low)
    hold on;
    plot(Posterior_all.marginals.normX*2*sqrt(2),Posterior_all.marginals.norm,'k')
    plot(Posterior_all.marginals.highX*2*sqrt(2),Posterior_all.marginals.high,'r');
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    xlabel('Distance from original position (cm)');
    ylabel('Decoded posterior');
    line(xlim, [0 0], 'linestyle','--','color','k');
    line([0 0], ylim, 'linestyle','--','color','k');
    
end
%% Calculating the field widths
Posterior_all.fieldWidths.norm = 2*(Posterior_all.decoder.numBins - calculate50crosses(Posterior_all.decoder.model.meanModel));
Posterior_all.fieldWidths.high = 2*(Posterior_all.decoder.numBins - calculate50crosses(Posterior_all.decoder_high.model.meanModel));
Posterior_all.fieldWidths.low =  2*(Posterior_all.decoder.numBins - calculate50crosses(Posterior_all.decoder_low.model.meanModel));

Posterior_all.fieldWidths.high_decPos = 2*(Posterior_all.decoder.numBins - calculate50crosses(Posterior_all.decoder_high_decPos.model.meanModel));
Posterior_all.fieldWidths.low_decPos =  2*(Posterior_all.decoder.numBins - calculate50crosses(Posterior_all.decoder_low_decPos.model.meanModel));
Posterior_all.fieldWidths.norm_decPos =  2*(Posterior_all.decoder.numBins - calculate50crosses(Posterior_all.decoder_norm_decPos.model.meanModel));
%%
if plot_all
    %     figure(15);
    %     set(15, 'Position', [100    100   1000   700]);
    f = figure;
    set(f, 'Position', [100    100   1000   700]);
    
    subplot(441)
    imagesc(Posterior_all.meanPost_new.low.early)
    text(3,-5,num2str(length(intersect(outcome.lowContrast, outcome.earlyTrials))));
    subplot(442)
    imagesc(Posterior_all.meanPost_new.low.correct)
    text(3,-5,num2str(length(intersect(outcome.lowContrast, outcome.correctTrials))));
    subplot(443)
    imagesc(Posterior_all.meanPost_new.low.late)
    text(3,-5,num2str(length(intersect(outcome.lowContrast, outcome.lateTrials))));
    subplot(444)
    imagesc(Posterior_all.meanPost_new.low.miss)
    text(3,-5,num2str(length(intersect(outcome.lowContrast, outcome.missTrials))));
    subplot(445)
    imagesc(Posterior_all.meanPost_new.norm.early)
    text(3,-5,num2str(length(intersect(outcome.normContrast, outcome.earlyTrials))));
    subplot(446)
    imagesc(Posterior_all.meanPost_new.norm.correct)
    text(3,-5,num2str(length(intersect(outcome.normContrast, outcome.correctTrials))));
    subplot(447)
    imagesc(Posterior_all.meanPost_new.norm.late)
    text(3,-5,num2str(length(intersect(outcome.normContrast, outcome.lateTrials))));
    subplot(448)
    imagesc(Posterior_all.meanPost_new.norm.miss)
    text(3,-5,num2str(length(intersect(outcome.normContrast, outcome.missTrials))));
    subplot(449)
    imagesc(Posterior_all.meanPost_new.high.early)
    text(3,-5,num2str(length(intersect(outcome.highContrast, outcome.earlyTrials))));
    subplot(4,4,10)
    imagesc(Posterior_all.meanPost_new.high.correct)
    text(3,-5,num2str(length(intersect(outcome.highContrast, outcome.correctTrials))));
    subplot(4,4,11)
    imagesc(Posterior_all.meanPost_new.high.late)
    text(3,-5,num2str(length(intersect(outcome.highContrast, outcome.lateTrials))));
    subplot(4,4,12)
    imagesc(Posterior_all.meanPost_new.high.miss)
    text(3,-5,num2str(length(intersect(outcome.highContrast, outcome.missTrials))));
    
    subplot(4,4,13)
    imagesc(Posterior_all.meanPost_new.all.early)
    text(3,-5,num2str(length(intersect(outcome.allContrast, outcome.earlyTrials))));
    subplot(4,4,14)
    imagesc(Posterior_all.meanPost_new.all.correct)
    text(3,-5,num2str(length(intersect(outcome.allContrast, outcome.correctTrials))));
    subplot(4,4,15)
    imagesc(Posterior_all.meanPost_new.all.late)
    text(3,-5,num2str(length(intersect(outcome.allContrast, outcome.lateTrials))));
    subplot(4,4,16)
    imagesc(Posterior_all.meanPost_new.all.miss)
    text(3,-5,num2str(length(intersect(outcome.allContrast, outcome.missTrials))));
    for n = 1:16
        subplot(4,4,n)
        axis xy; axis equal; axis tight;
        line([1 50], [1 50],'linestyle','--','color','k')
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        % set(gca,'CLim',[-0.5 0.5])
        set(gca, 'XTick',[0 25 50],'XTickLabel',[],'YTick',[0 25 50],'YTickLabel',[]);
    end
    subplot(4,4,1)
    ylabel('Low Contrast')
    title('EARLY')
    subplot(4,4,2)
    title('CORRECT')
    subplot(4,4,3)
    title('LATE')
    subplot(4,4,4)
    title('MISS')
    subplot(4,4,5)
    ylabel('Baseline Contrast');
    subplot(4,4,9)
    ylabel('High Contrast');
    subplot(4,4,13)
    ylabel('All');
    RedWhiteBlue;
    drawnow;
    pause(1);
end