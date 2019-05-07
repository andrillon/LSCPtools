%% initialisation

%% % parameters for reconstruction
fltrate=0.01;%sampled at 10Hz
timewindow=10;%size in seconds of the sliding window
stepslid=4;%step in seconds of the sliding window
maxdurationtrial=60;%duration maximum for trials for the sliding window (to be able to have the same size for each trial)

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

% try
%     load([ExpePath,filesep,Expe 'artifacts']);
%     channels_to_reject=define_channel_artifacts(artifact_var);%
% catch
    channels_to_reject=cell(1,length(subject_id));
% end

%get REM scored by hand
getpREMperiods

%%
for nS=1:length(subject_id)
    if isempty(pREM(nS).subID)==1
        continue
    end
    
    %%% Import in SPM
    subName=[Expe subject_id{nS}];
    SubID=subject_id{nS};
    fprintf('... Subject: %s\n',SubID);
    
    %%% Retrive behavioral data
    SubjectBehavData=load([BehavPath filesep subName filesep 'Result_' subName '.mat']);
    fprintf('... ... behavioral data retrived\n')
    
    %%% Retrieve events EEG data
    try
        load([ExtractedDataPath,'data_extracted' SubID '.mat'],'myevents')
    catch
        myevents=ft_read_event([rawEEGPath,subName,'.raw']);
    end
    
    %% retrieve structure of data : ordered by type of events, with
    saveName=[trialstPath,'trials_str',num2str(SubID)];
    CPR_trialevents
    CPR_trialstructure
    trial_sleep_scoring
    save(saveName,'trials_st')
    clear g

    %% training the model
    fprintf('... ... training model on di-otic\n')
    train_data=[];train_stim=[];
    trials_training = [trials_st(1:2).id];
    
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
    for nLag=1:length(lagvec)
        [g(:,nLag),rstim] = StimuliReconstruction (train_stim(:,1)', train_data, [], [], -lagvec(nLag));
    end
    
    fprintf('... ... test model on di-chotic\n')
    nSta=1;
    
    trials_phase=trials_st(4);
    
    %% Sleeptest trials
    
    for nT=1:length([trials_phase.id])
        
        nSta=trials_st(4).tagSleepScoring(nT);
        if isempty(nSta) || nSta<5
            continue
        end
        
%         phasicREMs=sum(trials_phase.trials(nT,1)<[pREM(nS).time]*500 & trials_phase.trials(nT,2)>[pREM(nS).time]*500);
        phasicREMs=trials_st(4).tagSleepScoring(nT)==6;
        count(nS)=count(nS)+1;
        nrREMtrial{nS}(count(nS))=phasicREMs;
        percphasic{nS}(count(nS))=nSta;
        
        nE=trials_phase.id(nT);
        nE_id{nS}(count(nS))=nE;
        %saveName for data
%         saveName=sprintf('%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2)); 
            saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        %load and reject artifacts if specified in channels_to_reject which
        %channels to reject
        filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
        
        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
            side_Tale{nS}(count(nS))=1;
        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
            side_Tale{nS}(count(nS))=2;
        end
        
        %reconstruction
        for nLag=1:length(lagvec)
            [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g(:,nLag), -lagvec(nLag));
            
            [test_rho1, ~]=corr(rstimdum',filtData.env(:,1));
            [test_rho2, ~]=corr(rstimdum',filtData.env(:,2));
            
            if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                tem_att{nS}(count(nS),nLag)=test_rho1;
                tem_ign{nS}(count(nS),nLag)=test_rho2;
            elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                tem_att{nS}(count(nS),nLag)=test_rho2;
                tem_ign{nS}(count(nS),nLag)=test_rho1;
            end
        end
        
        %half
        if ~isempty(begS_half)
            for nSeg=1:length(begS_half)
                for nLag=1:length(lagvec)
                    thisSeg=begS_half(nSeg):min([endS_half(nSeg),length(filtData.eeg)]);
                    [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg(thisSeg,:), g(:,nLag), -lagvec(nLag));
                    [test_rho1, ~]=corr(rstimdum',filtData.env(thisSeg,1));
                    [test_rho2, ~]=corr(rstimdum',filtData.env(thisSeg,2));
                    
                    if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                        tem_att_half{nS}(count(nS),nLag,nSeg)=test_rho1;
                        tem_ign_half{nS}(count(nS),nLag,nSeg)=test_rho2;
                    elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                        tem_att_half{nS}(count(nS),nLag,nSeg)=test_rho2;
                        tem_ign_half{nS}(count(nS),nLag,nSeg)=test_rho1;
                    end
                end
            end
        end
        
        %continuous
        if ~isempty(begS)
            for nSeg=1:length(begS)
                for nLag=1:length(lagvec)
                    thisSeg=begS(nSeg):min([endS(nSeg),length(filtData.eeg)]);
                    [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg(thisSeg,:), g(:,nLag), -lagvec(nLag));
                    [test_rho1, ~]=corr(rstimdum',filtData.env(thisSeg,1));
                    [test_rho2, ~]=corr(rstimdum',filtData.env(thisSeg,2));
                    
                    if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                        tem_att_Seg{nS}(count(nS),nLag,nSeg)=test_rho1;
                        tem_ign_Seg{nS}(count(nS),nLag,nSeg)=test_rho2;
                    elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                        tem_att_Seg{nS}(count(nS),nLag,nSeg)=test_rho2;
                        tem_ign_Seg{nS}(count(nS),nLag,nSeg)=test_rho1;
                    end
                end
            end
        end
    end
    
end