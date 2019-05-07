clear all
CPR_path_define_Thomas;

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
lagVec2=-(0:5:50);
freqVec=0.5:0.5:7;
% nSc=0;
for nS=1:length(subject_id)
    
    %%% Subj id
    subName=[Expe subject_id{nS}];
    SubID=subject_id{nS};
    fprintf('... Subject: %s\n',SubID);
    
    %%% Retrieve behavioral data
    SubjectBehavData=load([BehavPath filesep subName filesep 'Result_' subName '.mat']);
    fprintf('... ... behavioral data retrived\n')
    
    %%% Retrieve events and EEG data
    saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,1,filterband(1),filterband(2));
    if exist([rawPath,'data_extracted' SubID '.mat'])~=0 && exist([DataPath filesep saveName '.mat'])~=0
        load([rawPath,'data_extracted' SubID '.mat'],'myevents')
    else
        continue; %myevents=ft_read_event([rawEEGPath,subName,'.raw']);
    end
    
    %%% retrieve structure of data : ordered by type of events, with
    %%% scoring in tagSleep
    CPR_trialevents
    CPR_trialstructure
    trial_sleep_scoring
    for nFreq=1:length(freqVec)
        for nLag=1:length(lagVec2)
            %%% training the model
            clear g
            fprintf('... ... training model on di-otic\n')
            train_data=[];train_stim=[];
            trials_training = [trials_st(1:2).id];
            
            for nE=trials_training
                
                %saveName for data
                saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
                %load and reject artifacts if specified in channels_to_reject which
                %channels to reject
                filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
                filtData2=filtData;
                for nr=1:65
                    filtData2.eeg(:,nr)=bandpass(filtData.eeg(:,nr),100,freqVec(nFreq),freqVec(nFreq)+1,3);
                end
                for nr=1:2
                    filtData2.env(:,nr)=bandpass(filtData.env(:,nr),100,freqVec(nFreq),freqVec(nFreq)+1,3);
                end
                train_data=[train_data ; filtData2.eeg];
                train_stim=[train_stim ; filtData2.env];
            end
            [g,rstim] = StimuliReconstruction (train_stim(:,1)', train_data, [], [], -0:-1:lagVec2(nLag));
            
            fprintf('... ... test model on di-chotic\n')
            nSta=1;
            %%% reconstruction Wake
            %Forced Wake Trial
            countTr=0;
            for nE=[trials_st(3).id]
                countTr=countTr+1;
                fprintf('%g ... %g ... %g\n',nFreq,nLag,nE)
                
                count(nS,nSta)=count(nS,nSta)+1;
                nE_id{nS,nSta}(count(nS,nSta))=nE;
                if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                    side_Tale{nS,nSta}(count(nS,nSta))=1;
                elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                    side_Tale{nS,nSta}(count(nS,nSta))=2;
                end
                
                %saveName for data
                saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
                %load and reject artifacts if specified in channels_to_reject which
                %channels to reject
                filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
                filtData2=filtData;
                for nr=1:65
                    filtData2.eeg(:,nr)=bandpass(filtData.eeg(:,nr),100,freqVec(nFreq),freqVec(nFreq)+1,3);
                end
                for nr=1:2
                    filtData2.env(:,nr)=bandpass(filtData.env(:,nr),100,freqVec(nFreq),freqVec(nFreq)+1,3);
                end
                
                fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
                
                %stim reconstruction
                [gdum,rstimdum] = StimuliReconstruction ([], [], filtData2.eeg, g, -0:-1:lagVec2(nLag));
                [test_rho1, ~]=corr(rstimdum',filtData2.env(:,1));
                [test_rho2, ~]=corr(rstimdum',filtData2.env(:,2));
                
                %take into account the side of the real speech
                if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                    rho_att(nS,countTr,nLag,nFreq)=test_rho1;
                    rho_ign(nS,countTr,nLag,nFreq)=test_rho2;
                elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                    rho_att(nS,countTr,nLag,nFreq)=test_rho2;
                    rho_ign(nS,countTr,nLag,nFreq)=test_rho1;
                end
                decoding_att(nS,countTr,nLag,nFreq)=rho_att(nS,countTr,nLag,nFreq)>rho_ign(nS,countTr,nLag,nFreq);
            end
        end
    end
end