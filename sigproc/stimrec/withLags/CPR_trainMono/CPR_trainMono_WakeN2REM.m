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

NameStages={'W','N2','REM'};
Numstages={[0,1],2,[5,6]};%scoring reference for each stage
maxdiff=[];
count=zeros(length(subject_id),length(Numstages));
tem_att=cell(length(subject_id),length(Numstages));
tem_ign=cell(length(subject_id),length(Numstages));
side_Tale=cell(length(subject_id),length(Numstages));
nE_id=cell(length(subject_id),length(Numstages));
tem_att_Seg=cell(length(subject_id),length(Numstages));
tem_ign_Seg=cell(length(subject_id),length(Numstages));

% try
%     load([ExpePath,filesep,Expe 'artifacts']);
%     channels_to_reject=define_channel_artifacts(artifact_var);%
% catch
    channels_to_reject=cell(1,length(subject_id));
% end
%% 
for nS=length(subject_id)

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
    
    %% training the model
    clear g
    fprintf('... ... training model on di-otic\n')
    train_data=[];train_stim=[];
    trials_training = [trials_st(1:2).id];
    
    for nE=trials_training
        
        %saveName for data
            saveName=sprintf('%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
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
    %% reconstruction Wake
    
    %Forced Wake Trials
    for nE=[trials_st(3).id]
        
        count(nS,nSta)=count(nS,nSta)+1;
        nE_id{nS,nSta}(count(nS,nSta))=nE;
        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
            side_Tale{nS,nSta}(count(nS,nSta))=1;
        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
            side_Tale{nS,nSta}(count(nS,nSta))=2;
        end
        
        %saveName for data
            saveName=sprintf('%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        %load and reject artifacts if specified in channels_to_reject which
        %channels to reject
        filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
        
        fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
        for nLag=1:length(lagvec)
            
            %stim reconstruction
            [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g(:,nLag), -lagvec(nLag));
            [test_rho1, ~]=corr(rstimdum',filtData.env(:,1));
            [test_rho2, ~]=corr(rstimdum',filtData.env(:,2));
            
            %take into account the side of the real speech
            if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                tem_att{nS,nSta}(count(nS,nSta),nLag)=test_rho1;
                tem_ign{nS,nSta}(count(nS,nSta),nLag)=test_rho2;
            elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                tem_att{nS,nSta}(count(nS,nSta),nLag)=test_rho2;
                tem_ign{nS,nSta}(count(nS,nSta),nLag)=test_rho1;
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
                        tem_att_half{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho1;
                        tem_ign_half{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho2;
                    elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                        tem_att_half{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho2;
                        tem_ign_half{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho1;
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
                        tem_att_Seg{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho1;
                        tem_ign_Seg{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho2;
                    elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                        tem_att_Seg{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho2;
                        tem_ign_Seg{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho1;
                    end
                end
            end
        end
    end
    
    %% reconstruction Sleep

    % Sleeptest trials
    for nT=1:length([trials_st(4).id])
        %select N2 and REM
        nSta=find(cellfun(@(x) (trials_st(4).tagSleepScoring(nT)>=x(1) && trials_st(4).tagSleepScoring(nT)<=x(end)),Numstages(2:end)))+1;
        if isempty(nSta)
            continue
        end
        count(nS,nSta)=count(nS,nSta)+1;
        nE=trials_st(4).id(nT);
        nE_id{nS,nSta}(count(nS,nSta))=nE;
        
            saveName=sprintf('%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
        
        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
            side_Tale{nS,nSta}(count(nS,nSta))=1;
        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
            side_Tale{nS,nSta}(count(nS,nSta))=2;
        end
        
        fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
        for nLag=1:length(lagvec)
            [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g(:,nLag), -lagvec(nLag));
            
            [test_rho1, ~]=corr(rstimdum',filtData.env(1:length(filtData.eeg),1));
            [test_rho2, ~]=corr(rstimdum',filtData.env(1:length(filtData.eeg),2));
            
            if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                tem_att{nS,nSta}(count(nS,nSta),nLag)=test_rho1;
                tem_ign{nS,nSta}(count(nS,nSta),nLag)=test_rho2;
            elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                tem_att{nS,nSta}(count(nS,nSta),nLag)=test_rho2;
                tem_ign{nS,nSta}(count(nS,nSta),nLag)=test_rho1;
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
                        tem_att_half{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho1;
                        tem_ign_half{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho2;
                    elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                        tem_att_half{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho2;
                        tem_ign_half{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho1;
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
                        tem_att_Seg{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho1;
                        tem_ign_Seg{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho2;
                    elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                        tem_att_Seg{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho2;
                        tem_ign_Seg{nS,nSta}(count(nS,nSta),nLag,nSeg)=test_rho1;
                    end
                end
            end
        end
        
    end
    save([savePath,filesep,'Reconstruction',num2str(filterband(1)),'to',num2str(filterband(2)),'.mat'],'tem_att','tem_att_half','tem_ign_half','tem_ign','tem_att_Seg','subject_id','tem_ign_Seg','begS','endS','begS_half','endS_half','subject_id','count','lagvec','Numstages','NameStages','side_Tale','nE_id','userID','Expe','channels_to_reject')%% plot for subjects
end