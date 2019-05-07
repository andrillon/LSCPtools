%%% set parameters for reconstruction
lagvec=[0:50];% lag for linear model
restrictREM='all';%'burst','isolated'

%get REM scored by hand
getpREMperiods
REMtonic_all=cell(1,length(subject_id));
REMphasic_all=cell(1,length(subject_id));

% try
%     load([ExpePath,filesep,Expe 'artifacts']);
%     channels_to_reject=define_channel_artifacts(artifact_var);%
% catch
    channels_to_reject=cell(1,length(subject_id));
% end
%%
for nS=1:length(subject_id)
    
    if isempty(pREM(nS).subID)==1
        continue
    end
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
    
    %% training the model
    trials_training = [trials_st(1:2).id];
    train_data=[];train_stim=[];
    for nE=trials_training
        %      fprintf('.. %g/40 ..\n',nE)
            saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        
        filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
        %         filtData.eeg=(filtData.eeg-repmat(min(filtData.eeg),[size(filtData.eeg,1) 1]))./repmat(max(filtData.eeg)-min(filtData.eeg),[size(filtData.eeg,1) 1]);
        
        train_data=[train_data;filtData.eeg(1:length(filtData.eeg),:)];
        train_stim=[train_stim;filtData.env(1:length(filtData.eeg),:)];
    end
    for nLag=1:length(lagvec)
        [g(:,nLag),rstim] = StimuliReconstruction(train_stim(:,1)', train_data, [], [], -lagvec(nLag));
    end
    
    %%
    ExpPhase=4;%select only sleeptest trials
    trials_phase=trials_st(ExpPhase);
    if strcmp(restrictREM,'burst')
        REMper=pREM(nS).periodsBurstREM;
    elseif strcmp(restrictREM,'isolated')
        REMper=pREM(nS).periodsIsolatedREM;
    elseif strcmp(restrictREM,'all')
        REMper=pREM(nS).periodsAllREM;
    end
    
    %%
    if isempty(REMper)
        continue
    end
    REMtonic_eeg=[];REMtonic_stim=[];
    REMphasic_eeg=[];REMphasic_stim=[];
    
    %pick all trials scored as REM
    trialsREM=find(trials_phase.tagSleepScoring>4);
    %%
    for nrtrial_REM=trialsREM
        
        %%%is rem within a trial ?
        trialrem=find((REMper(:,1)>=trials_phase.trials(nrtrial_REM,1)/500 & REMper(:,1)<trials_phase.trials(nrtrial_REM,2)/500)...
            | (REMper(:,2)>=trials_phase.trials(nrtrial_REM,1)/500 & REMper(1:end,2)<trials_phase.trials(nrtrial_REM,2)/500));
        %%
        %             try
        nE=trials_phase.id(nrtrial_REM);
        clear filtData
            saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
        eegdata=filtData.eeg;
        
        %%%get REM periods within the trial
        REMphasic_trial=[];
        for rem=trialrem'
            startrem=floor(max([1,[REMper(rem,1)-trials_phase.trials(nrtrial_REM,1)/500]*100]));
            endrem=floor(min([length(eegdata),[REMper(rem,2)-trials_phase.trials(nrtrial_REM,1)/500]*100]));
            REMphasic_trial=[REMphasic_trial,startrem:endrem];
        end
        
        %%%concatenate periods of pREM and tREM within a trial scored as
        %%%REM
        REMtonic_trial = setdiff(1:length(eegdata),REMphasic_trial);
        REMtonic_all{nS} = [REMtonic_all{nS},REMtonic_trial];
        REMphasic_all{nS} = [REMphasic_all{nS},REMphasic_trial];
        REMtonic_eeg = [REMtonic_eeg;filtData.eeg(REMtonic_trial,:)];
        REMphasic_eeg = [REMphasic_eeg;filtData.eeg(REMphasic_trial,:)];
        
        if SubjectBehavData.TrialsCaracs(trials_phase.id(nrtrial_REM)).tSide=='L'
            REMtonic_stim = [REMtonic_stim;filtData.env(REMtonic_trial,:)];
            REMphasic_stim = [REMphasic_stim;filtData.env(REMphasic_trial,:)];
        elseif SubjectBehavData.TrialsCaracs(trials_phase.id(nrtrial_REM)).tSide=='R'
            REMtonic_stim = [REMtonic_stim;fliplr(filtData.env(REMtonic_trial,:))];
            REMphasic_stim = [REMphasic_stim;fliplr(filtData.env(REMphasic_trial,:))];
        end
    end
    
    %%
    lenphases(nS,1)=size(REMtonic_stim,1);
    lenphases(nS,2)=size(REMphasic_stim,1);
    for nLag=1:length(lagvec)
        if size(REMtonic_stim,1)>0
            [gdum,rstim_tonic] = StimuliReconstruction ([], [], REMtonic_eeg, g(:,nLag), -lagvec(nLag));
            [rho_real_tonic(nS,nLag), ~]=corr(rstim_tonic',REMtonic_stim(:,1));
            [rho_jab_tonic(nS,nLag), ~]=corr(rstim_tonic',REMtonic_stim(:,2));
        end
        if size(REMphasic_stim,1)>0
            [gdum,rstim_phasic] = StimuliReconstruction ([], [], REMphasic_eeg, g(:,nLag), -lagvec(nLag));
            [rho_real_phasic(nS,nLag), ~]=corr(rstim_phasic',REMphasic_stim(:,1));
            [rho_jab_phasic(nS,nLag), ~]=corr(rstim_phasic',REMphasic_stim(:,2));
        end
    end
    
end
