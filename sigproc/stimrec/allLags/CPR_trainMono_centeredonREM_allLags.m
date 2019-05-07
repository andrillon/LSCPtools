%% Parameters for sliding window
fltrate=0.01;
smallEpoch = 4;               % Size of the window for the computation of pearson's R between reconstruction and real streams (correlation window)
flankSize = 10;               % Data used on both sides of the micro-event
stepslid=0.1;                   %step for sliding window
fs=500;%Hz
getpREMperiods

ncountREM=zeros(1,length(subject_id));
REM_valid=cell(1,length(subject_id));
REM_trial=cell(1,length(subject_id));
timewindow = smallEpoch;%sliding window ? : if empty, entire trial is taken
maxdurationtrial=flankSize*2.01+smallEpoch;%to compare btw trials
%function to get the segments beginning and ends and lengths for the
%sliding window
[begS,endS,lenS]=time_window_parameter(timewindow,stepslid,maxdurationtrial,fltrate);
channels_to_reject=cell(1,length(subject_id));

%%

for nS=1:length(subject_id)
    
    
     %%% Subj id
    subName=[Expe subject_id{nS}];
    SubID=subject_id{nS};
    fprintf('... Subject: %s\n',SubID);
    
    %%% Retrieve behavioral data
    SubjectBehavData=load([BehavPath filesep subName filesep 'Result_' subName '.mat']);
    fprintf('... ... behavioral data retrived\n')
    
    %%% Retrieve events and EEG data
    try
        load([ExtractedDataPath,'data_extracted' SubID '.mat'],'myevents')
    catch
        myevents=ft_read_event([rawEEGPath,subName,'.raw']);
    end    
    
    %%
    CPR_trialevents
    CPR_trialstructure
    trial_sleep_scoring
    clear g
    
    %% training
    trials_training = [trials_st(1:2).id];
    train_data=[];train_stim=[];

    try
        load([model_path subName,num2str(filterband(1)),'to',num2str(filterband(2)),'.mat'],'g')
    catch
        for nE=trials_training
            
            %saveName for data
            %         saveName=sprintf('%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
            saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
            %load and reject artifacts if specified in channels_to_reject which
            %channels to reject
            filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
            
            train_data=[train_data ; filtData.eeg];
            train_stim=[train_stim ; filtData.env];
        end
        [g,rstim] = StimuliReconstruction (train_stim(:,1)', train_data, [], [], lagvec);
        save([model_path subName,num2str(filterband(1)),'to',num2str(filterband(2)),'.mat'],'g')
    end    
    %%
    ExpPhase=4;
    trials_phase = trials_st(ExpPhase);
    REMscaracs=pREM(nS).caracs;
    if isempty(REMscaracs)
        continue
    end
    
    REMstart=REMscaracs.start;
    REMtime=pREM(nS).time;
    ncountREM(nS)=0;
    %%
    for rem=1:length(REMstart)
        %             try
        rem_selec = REMstart(rem);
        remtime = REMtime(rem_selec)*fs;%Hz
        %is rem within a trial ?
        
        trialrem=remtime-(flankSize+smallEpoch*2.1)*fs>trials_phase.trials(:,1) & remtime+ceil((flankSize+smallEpoch*2.1)*fs)<trials_phase.trials(:,2);
        keyboard
        if sum(trialrem)==0
            continue
        else
            %%
            %             try
            if trials_phase.tagSleepScoring(trialrem)<5;
                continue
            end
            clear filtData
            nrtrial_REM=find(trialrem);
            nE=trials_phase.id(trialrem);
            %pick all trials scored as REM
            
            lockvalue=[remtime-trials_phase.trials(nrtrial_REM,1)-ceil((flankSize+smallEpoch/2)*fs)]/(fs*fltrate);

            saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
            %load and reject artifacts if specified in channels_to_reject which
            %channels to reject
            filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
            eegdata=filtData.eeg;
            if floor(lockvalue+endS(end))>length(eegdata)
                continue
            end
            ncountREM(nS)=ncountREM(nS)+1;
            REM_valid{nS}(ncountREM(nS))=rem_selec;
            REM_trial{nS}(ncountREM(nS))=nrtrial_REM;
            
            SidenE=SubjectBehavData.TrialsCaracs(nE).tSide;
            for nSeg=1:length(begS)
                thisSeg=floor(lockvalue+[begS(nSeg):endS(nSeg)]);
                
                    [gdum,rstim] = StimuliReconstruction ([], [], eegdata(thisSeg,:), g, lagvec);
                    [test_rho1, ~]=corr(rstim',filtData.env(thisSeg,1));
                    [test_rho2, ~]=corr(rstim',filtData.env(thisSeg,2));
                    if SidenE=='L' %R
                        %                         side_Tale_Sleep(nS,nEc)=1;
                        tem_att{nS}(ncountREM(nS),nSeg)=test_rho1;
                        tem_ign{nS}(ncountREM(nS),nSeg)=test_rho2;
                    elseif SidenE=='R' %L
                        %                         side_Tale_Sleep(nS,nEc)=2;
                        tem_att{nS}(ncountREM(nS),nSeg)=test_rho2;
                        tem_ign{nS}(ncountREM(nS),nSeg)=test_rho1;
                    end
            end
        end
    end
    %trial nr from sleep new trial,remnr for this subject,start time,start time,burst?,left?,right?,up?,down?
    REM_compute{nS}=[REM_trial{nS}(:),REM_valid{nS}(:),REMtime(REM_valid{nS})',REMtime(REM_valid{nS})',ismember(REM_valid{nS},REMscaracs.burst)',ismember(REM_valid{nS}-1,[0,REMscaracs.end])',ismember(REM_valid{nS},REMscaracs.left)',ismember(REM_valid{nS},REMscaracs.right)',ismember(REM_valid{nS},REMscaracs.up)',ismember(REM_valid{nS},REMscaracs.down)'];
    save([savePath,'reconstruction_REMcentered',num2str(filterband(1)),'to',num2str(filterband(2)),'.mat'],'REM_compute','subject_id','tem_att','tem_ign','REM_trial','REM_valid','REM_compute','begS','endS','smallEpoch','flankSize','lagvec','userID','Expe','channels_to_reject')    
end