function [es, spikeRate, VRdata, VRdata_o, chans, spontRate, spikeTimes] = ...
    getVRspikes(animal,iseries,iexp, ...
    discreet_steps, flag_load, flag_spont, flag_spkrate, igroup,iaddinfo,loadTimes)

%% initialization
global DIRS
global THRES

SetDefaultDirs

if nargin~=5
    flag_load = 1;
end

if nargin<4
    discreet_steps = 100;
end

if nargin<6
    flag_spont = 1;
end

if nargin<7
    flag_spkrate = 0;
end

if nargin<10
    loadTimes = 0;
end

fname = [animal '_' num2str(iseries) '_' num2str(iexp)];
dDIRname = [DIRS.multichanspikes filesep animal filesep num2str(iseries)];

if flag_load
    if exist([dDIRname filesep fname '_contspikes' '.mat'],'file');
        try
            load([dDIRname filesep fname '_contspikes','.mat'],'es');
%         if exist('es','var')
            if size(es.trajpos,2)==discreet_steps
                run_rest = 0;
            else
                display('Old structure has different position sampling')
                run_rest = 1;
            end
        catch %else
            display('No even sampling')
            run_rest = 1;
        end
    else
        display('No existing file');
        run_rest = 1;
    end
else
    run_rest = 1;
end

if run_rest
    
    load([dDIRname filesep fname '_screenTimes']);
    
    [VRdata, VRdata_o, es_test] = VRWheelLoad(animal, iseries, iexp);
    try
        if nargin<8
            inp = inputdlg('Enter the group to be loaded:','0');
            igroup = str2num(inp{1});            
        end
        for igroup_idx = 1:length(igroup)
            if exist('iaddinfo')
                temp = getKwikSpikes(animal, iseries, iexp, igroup(igroup_idx),iaddinfo);
            else
                temp = getKwikSpikes(animal, iseries, iexp, igroup(igroup_idx));
            end
           tempNumCells = length(temp);
           if igroup_idx==1
               chans=temp;
               numCells = tempNumCells;
           else
               for icell = 1:tempNumCells
                   chans(numCells+1) = temp(icell);
                   numCells = numCells + 1;
               end
           end
        end
        isolDist = zeros(1,length(chans));
    catch
        chans = UnitLoad([DIRS.spikes filesep 'Klustered'], animal, iseries, iexp);
        isolDist = zeros(1,length(chans));
        for icellIsol = 1:length(chans)
            if ~isempty(chans(icellIsol).isolDist)
                isolDist(icellIsol) = chans(icellIsol).isolDist(1);
            else
                isolDist(icellIsol) = 0;
            end
        end
    end
    
    screenTimes = screenTimes./30000;
    %  Correcting and matching the recording and the VR
    RecToVR_correction = screenTimes(1) - es_test.screenTimes2(2);
    screenTimes = screenTimes - RecToVR_correction;
    
    %% To get the even sampled parameters
    
    evenSampleTime = (1/60):(1/60):max(screenTimes);
    
    es = VREvenSample(VRdata, evenSampleTime, screenTimes);
    es.iexp = zeros(size(es.sampleTimes));
    es.iexp(:) = iexp;
    
    es.isolDist = isolDist;
    
    % getting the spike trains
    numCells = length(chans);
    es.spikeTrain = zeros(length(es.sampleTimes), numCells);
    for ichan = 1:length(chans)
        chans(ichan).spiketimes = chans(ichan).spiketimes - RecToVR_correction;
        chans(ichan).spiketimes(chans(ichan).spiketimes<0) = [];
        
%         display(['Processing cell: ' num2str(ichan)]);
        % finding the bin for each spike
        st = ceil(chans(ichan).spiketimes./(1/60));
        numExcessSpikes = sum(st>length(es.sampleTimes));
        st(st>length(es.sampleTimes)) = [];
        st(st==0) = [];
        
        if numExcessSpikes>0
            chans(ichan).spiketimes(end-numExcessSpikes+1:end) = [];
        end
        es.spikeTrain(st,ichan) = 1;
        
        repeats = st(diff(st)==0);
        for irep = 1:length(repeats)
            es.spikeTrain(repeats(irep),ichan) = es.spikeTrain(repeats(irep),ichan) + 1;
        end
        es.spikeIDs{ichan} = chans(ichan).id;
        if loadTimes
            es.spikeTimes{ichan} = chans(ichan).spiketimes;
        end
        clear st trash st_unique order repeats irep
    end
    es.mua = sum(es.spikeTrain,2);
    
%     try
%         traj    = VRdata.TRIAL.traj(1:(max(find(VRdata.TRIAL.traj))-1+VRdata.EXP.pause_frames));
%     catch
%         % this is for when the experiment is quit in between a run
%         display('This experiment was quit between a run!');
%         traj    = VRdata.TRIAL.traj(1:(max(find(VRdata.TRIAL.traj))-1));
%     end
    
    %% Checking the the sync between the display and recording
    if length(es.screenTimes) ~= length(es.screenTimes2)-1
        disp(['!!!!!WARNING!!!!!!: ' num2str(animal) '_' num2str(iseries) '_' num2str(iexp) ' trajectory has ' num2str(length(es.screenTimes) - length(es.screenTimes2)+1) ' more frames than screen refreshes']);
    end
    
    
    %% get all spikeTimes
    if flag_spkrate
        for ichan = 1:length(chans)
            if isempty(chans(ichan).spiketimes)
                spikeTimes(ichan).t = [];
                spikeRate(ichan).t = [];
            else
                spikeTimes(ichan).t = chans(ichan).spiketimes;
                for iscreen = 2:length(screenTimes)
                    nspike = sum((spikeTimes(ichan).t > screenTimes(iscreen-1)) & (spikeTimes(ichan).t <= screenTimes(iscreen)));
                    spikeRate(ichan).t(iscreen-1) = nspike/(screenTimes(iscreen)-screenTimes(iscreen-1));
                end
            end
            
        end
    end
    %% the spont rate for the period before the start of the VR
    if flag_spont
        spontRate = zeros(1,length(chans));
        for cellID = 1:length(chans)
            numSpks(cellID) = sum(chans(cellID).spiketimes<screenTimes(1));
            
        end
        spontRate = numSpks./screenTimes(1);
        clear numSpks
    else
        spontRate = [];
    end
    
    %% save them to expt
    if exist([dDIRname filesep fname '_contspikes' '.mat'],'file');
        display('Appending the new structures to the file!');
        if flag_spkrate
            save([dDIRname filesep fname '_contspikes'],'spikeRate','spikeTimes','es', '-append');
        else
            save([dDIRname filesep fname '_contspikes'],'es', '-append');
        end
    else
        display('Saving the file!');
        display('Problem with append, saving without append!')
        
        if flag_spkrate
            save([dDIRname filesep fname '_contspikes'],'spikeRate','VRdata','VRdata_o','chans','spikeTimes','es');
        else
            save([dDIRname filesep fname '_contspikes'],'VRdata','VRdata_o','chans','es');
        end
    end
    save([dDIRname filesep fname '_es'],'es');
end

