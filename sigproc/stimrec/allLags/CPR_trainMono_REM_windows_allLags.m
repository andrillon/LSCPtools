%function to get the segments beginning and ends and lengths for the
%sliding window
% [begS,endS,lenS]=time_window_parameter(timewindow,stepslid,maxdurationtrial,fltrate);
begS=[];endS=[];
%time windows for 0-30s,30s-60s
begS_half=[1,30./fltrate];endS_half=[30/fltrate+1,60/fltrate];

%initialization
maxdiff=[];
count=zeros(1,length(subject_id));
tem_att=cell(1,length(subject_id));
tem_ign=cell(1,length(subject_id));
tem_att_half=cell(1,length(subject_id));
tem_ign_half=cell(1,length(subject_id));
tem_att_Seg=cell(1,length(subject_id));
tem_ign_Seg=cell(1,length(subject_id));
side_Tale=cell(1,length(subject_id));
nrREMtrial=cell(1,length(subject_id));
percphasic=cell(1,length(subject_id));
nE_id=cell(1,length(subject_id));

channels_to_reject=cell(1,length(subject_id));

for nS=1:length(subject_id)
    
    SubID=subject_id{nS};
    subName=[Expe SubID];
    
    %%% get behavioral data
    SubjectBehavData = load([BehavPath filesep subName filesep 'Result_' subName '.mat']);
    
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
    clear g
    
    if sum(hypnogram>4)==0
        continue
    end
    
    %% training the model
    trials_training = [trials_st(1:2).id];

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
    ExpPhase=4;%select only sleeptest trials
    trials_phase=trials_st(ExpPhase);
    trialsREM=find(trials_phase.tagSleepScoring>4);
    
    %%
    for nT=1:length([trials_phase.id])
        
        nSta=trials_st(4).tagSleepScoring(nT);
        if isempty(nSta) || nSta<5
            continue
        end
        
        nE=trials_phase.id(nT);
        
        %find windows of scoring within the trial
        windowsREM=find(time_scoring(1:end-1)>=trials_phase.trials(nT,1)/500 & time_scoring(2:end)<trials_phase.trials(nT,2)/500);
        clear filtData
        
        saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
        eegdata=filtData.eeg;
        
        
        for tr=windowsREM
            
            phaserem=hypnogram(tr);
            count(nS)=count(nS)+1;
            phaseREM{nS}(count(nS))=phaserem;
            nE_id{nS}(count(nS))=nE;
            
            %%%get period of window
            startwindow=floor((time_scoring(tr)-trials_phase.trials(nT,1)/500)*100);
            endwindow=floor((time_scoring(tr+1)-trials_phase.trials(nT,1)/500)*100);
            
            if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                side_Tale{nS}(count(nS))=1;
            elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                side_Tale{nS}(count(nS))=2;
            end
            
            %reconstruction
            [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg(startwindow:endwindow,:), g, lagvec);
            
            [test_rho1, ~]=corr(rstimdum',filtData.env(startwindow:endwindow,1));
            [test_rho2, ~]=corr(rstimdum',filtData.env(startwindow:endwindow,2));
            
            if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                tem_att{nS}(count(nS))=test_rho1;
                tem_ign{nS}(count(nS))=test_rho2;
            elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                tem_att{nS}(count(nS))=test_rho2;
                tem_ign{nS}(count(nS))=test_rho1;
            end
        end
    end
end