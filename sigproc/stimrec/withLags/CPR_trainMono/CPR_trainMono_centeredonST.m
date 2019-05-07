%% Parameters for sliding window
fltrate=0.01;
smallEpoch = 4;               % Size of the window for the computation of pearson's R between reconstruction and real streams (correlation window)
flankSize = 10;               % Data used on both sides of the micro-event
stepslid=0.1;                   %step for sliding window
fs=500;%Hz

ncountST=zeros(1,length(subject_id));
ST_valid=cell(1,length(subject_id));
ST_trial=cell(1,length(subject_id));
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
    
    %%% retrieve structure of data : ordered by type of events, with
    %%% scoring in tagSleep
    CPR_trialevents
    CPR_trialstructure
    trial_sleep_scoring
    
    
    %%
    saveName=[trialstPath,'trials_str',num2str(SubID)];
    CPR_trialevents
    CPR_trialstructure
    trial_sleep_scoring
    clear g
    
    %% training
    trials_training = [trials_st(1:2).id];
    train_data=[];train_stim=[];
    for nE=trials_training
        %      fprintf('.. %g/40 ..\n',nE)
        %         saveName=sprintf('%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        %load and reject artifacts if specified in channels_to_reject which
        %channels to reject
        filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
        
        train_data=[train_data;filtData.eeg];
        train_stim=[train_stim;filtData.env];
    end
    for nLag=1:length(lagvec)
        [g(:,nLag),rstim] = StimuliReconstruction(train_stim(:,1)', train_data, [], [], -lagvec(nLag));
    end
    
    %%
    ExpPhase=4;
    trials_phase = trials_st(ExpPhase);
    load([savePath, 'ST_',subName,'.mat']);
    if isempty(ST)
        continue
    end
    STnegpeak=[ST.CWT_NegativePeak];
    STtime=STnegpeak/Info.Recording.sRate;
    ncountST(nS)=0;
    %%
    for st_idx=1:length(STnegpeak)
        %             try
        rem_selec = st_idx;
        st_time = STnegpeak(st_idx);%Hz
        %is st_idx within a trial ?
        
        trialst=st_time-(flankSize+smallEpoch*2.1)*fs>trials_phase.trials(:,1) & st_time+ceil((flankSize+smallEpoch*2.1)*fs)<trials_phase.trials(:,2);
        
        if sum(trialst)==0
            continue
        else
            %%
            %             try
            if trials_phase.tagSleepScoring(trialst)<5;
                continue
            end
            clear filtData
            nrtrial_ST=find(trialst);
            nE=trials_phase.id(trialst);
            %pick all trials scored as ST
            
            lockvalue=[st_time-trials_phase.trials(nrtrial_ST,1)-ceil((flankSize+smallEpoch/2)*fs)]/(fs*fltrate);

            saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
            %load and reject artifacts if specified in channels_to_reject which
            %channels to reject
            filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
            eegdata=filtData.eeg;
            if floor(lockvalue+endS(end))>length(eegdata)
                continue
            end
            ncountST(nS)=ncountST(nS)+1;
            ST_valid{nS}(ncountST(nS))=rem_selec;
            ST_trial{nS}(ncountST(nS))=nrtrial_ST;
            
            SidenE=SubjectBehavData.TrialsCaracs(nE).tSide;
            for nSeg=1:length(begS)
                thisSeg=floor(lockvalue+[begS(nSeg):endS(nSeg)]);
                
                for nLag=1:length(lagvec)
                    [gdum,rstim] = StimuliReconstruction ([], [], eegdata(thisSeg,:), g(:,nLag), -lagvec(nLag));
                    [test_rho1, ~]=corr(rstim',filtData.env(thisSeg,1));
                    [test_rho2, ~]=corr(rstim',filtData.env(thisSeg,2));
                    if SidenE=='L' %R
                        %                         side_Tale_Sleep(nS,nEc)=1;
                        tem_att{nS}(ncountST(nS),nLag,nSeg)=test_rho1;
                        tem_ign{nS}(ncountST(nS),nLag,nSeg)=test_rho2;
                    elseif SidenE=='R' %L
                        %                         side_Tale_Sleep(nS,nEc)=2;
                        tem_att{nS}(ncountST(nS),nLag,nSeg)=test_rho2;
                        tem_ign{nS}(ncountST(nS),nLag,nSeg)=test_rho1;
                    end
                end
            end
        end
    end
    %trial nr from sleep new trial,remnr for this subject,start time,start time,burst?,left?,right?,up?,down?
    ST_compute{nS}=[ST_trial{nS}(:),ST_valid{nS}(:),STtime(ST_valid{nS})',STtime(ST_valid{nS})',ismember(ST_valid{nS},STscaracs.burst)',ismember(ST_valid{nS}-1,[0,STscaracs.end])',ismember(ST_valid{nS},STscaracs.left)',ismember(ST_valid{nS},STscaracs.right)',ismember(ST_valid{nS},STscaracs.up)',ismember(ST_valid{nS},STscaracs.down)'];
    save([savePath,'reconstruction_STcentered',num2str(filterband(1)),'to',num2str(filterband(2)),'.mat'],'ST_compute','subject_id','tem_att','tem_ign','ST_trial','ST_valid','ST_compute','begS','endS','smallEpoch','flankSize','lagvec','userID','Expe','channels_to_reject')    
end