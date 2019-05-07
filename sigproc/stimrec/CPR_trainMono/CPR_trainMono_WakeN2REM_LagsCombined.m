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
lagVec2=-(1:50);
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
    
    %% training the model
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
        
        train_data=[train_data ; filtData.eeg];
        train_stim=[train_stim ; filtData.env];
    end
    [g,rstim] = StimuliReconstruction (train_stim(:,1)', train_data, [], [], lagVec2);
    
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
        saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        %load and reject artifacts if specified in channels_to_reject which
        %channels to reject
        filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
        
        fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
        
        %stim reconstruction
        [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g, lagVec2);
        [test_rho1, ~]=corr(rstimdum',filtData.env(:,1));
        [test_rho2, ~]=corr(rstimdum',filtData.env(:,2));
        
        %take into account the side of the real speech
        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
            tem_att{nS,nSta}(count(nS,nSta))=test_rho1;
            tem_ign{nS,nSta}(count(nS,nSta))=test_rho2;
        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
            tem_att{nS,nSta}(count(nS,nSta))=test_rho2;
            tem_ign{nS,nSta}(count(nS,nSta))=test_rho1;
        end
        
        %half
        if ~isempty(begS_half)
            for nSeg=1:length(begS_half)
                thisSeg=begS_half(nSeg):min([endS_half(nSeg),length(filtData.eeg)]);
                [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg(thisSeg,:), g, lagVec2);
                [test_rho1, ~]=corr(rstimdum',filtData.env(thisSeg,1));
                [test_rho2, ~]=corr(rstimdum',filtData.env(thisSeg,2));
                
                if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                    tem_att_half{nS,nSta}(count(nS,nSta),nSeg)=test_rho1;
                    tem_ign_half{nS,nSta}(count(nS,nSta),nSeg)=test_rho2;
                elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                    tem_att_half{nS,nSta}(count(nS,nSta),nSeg)=test_rho2;
                    tem_ign_half{nS,nSta}(count(nS,nSta),nSeg)=test_rho1;
                end
            end
        end
        
        %continuous
        if ~isempty(begS)
            for nSeg=1:length(begS)
                thisSeg=begS(nSeg):min([endS(nSeg),length(filtData.eeg)]);
                [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg(thisSeg,:), g, lagVec2);
                [test_rho1, ~]=corr(rstimdum',filtData.env(thisSeg,1));
                [test_rho2, ~]=corr(rstimdum',filtData.env(thisSeg,2));
                
                if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                    tem_att_Seg{nS,nSta}(count(nS,nSta),nSeg)=test_rho1;
                    tem_ign_Seg{nS,nSta}(count(nS,nSta),nSeg)=test_rho2;
                elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                    tem_att_Seg{nS,nSta}(count(nS,nSta),nSeg)=test_rho2;
                    tem_ign_Seg{nS,nSta}(count(nS,nSta),nSeg)=test_rho1;
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
        
        saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
        
        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
            side_Tale{nS,nSta}(count(nS,nSta))=1;
        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
            side_Tale{nS,nSta}(count(nS,nSta))=2;
        end
        
        fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
        [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g, lagVec2);
        
        [test_rho1, ~]=corr(rstimdum',filtData.env(1:length(filtData.eeg),1));
        [test_rho2, ~]=corr(rstimdum',filtData.env(1:length(filtData.eeg),2));
        
        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
            tem_att{nS,nSta}(count(nS,nSta))=test_rho1;
            tem_ign{nS,nSta}(count(nS,nSta))=test_rho2;
        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
            tem_att{nS,nSta}(count(nS,nSta))=test_rho2;
            tem_ign{nS,nSta}(count(nS,nSta))=test_rho1;
        end
        
        %half
        if ~isempty(begS_half)
            for nSeg=1:length(begS_half)
                thisSeg=begS_half(nSeg):min([endS_half(nSeg),length(filtData.eeg)]);
                [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg(thisSeg,:), g, lagVec2);
                [test_rho1, ~]=corr(rstimdum',filtData.env(thisSeg,1));
                [test_rho2, ~]=corr(rstimdum',filtData.env(thisSeg,2));
                
                if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                    tem_att_half{nS,nSta}(count(nS,nSta),nSeg)=test_rho1;
                    tem_ign_half{nS,nSta}(count(nS,nSta),nSeg)=test_rho2;
                elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                    tem_att_half{nS,nSta}(count(nS,nSta),nSeg)=test_rho2;
                    tem_ign_half{nS,nSta}(count(nS,nSta),nSeg)=test_rho1;
                end
            end
        end
        
        %continuous
        if ~isempty(begS)
            for nSeg=1:length(begS)
                thisSeg=begS(nSeg):min([endS(nSeg),length(filtData.eeg)]);
                [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg(thisSeg,:), g, lagVec2);
                [test_rho1, ~]=corr(rstimdum',filtData.env(thisSeg,1));
                [test_rho2, ~]=corr(rstimdum',filtData.env(thisSeg,2));
                
                if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                    tem_att_Seg{nS,nSta}(count(nS,nSta),nSeg)=test_rho1;
                    tem_ign_Seg{nS,nSta}(count(nS,nSta),nSeg)=test_rho2;
                elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                    tem_att_Seg{nS,nSta}(count(nS,nSta),nSeg)=test_rho2;
                    tem_ign_Seg{nS,nSta}(count(nS,nSta),nSeg)=test_rho1;
                end
            end
        end
        
    end
    nSc=nSc+1;
    for nSta=1:3
        mean_att(nSc,nSta)=mean(tem_att{nS,nSta});
        mean_ign(nSc,nSta)=mean(tem_ign{nS,nSta});
        mean_decoding(nSc,nSta)=mean(tem_att{nS,nSta}>tem_ign{nS,nSta});
    end
    save([savePath,filesep,'Reconstruction',num2str(filterband(1)),'to',num2str(filterband(2)),'_TA.mat'],'tem_att','tem_att_half','tem_ign_half','tem_ign','tem_att_Seg','subject_id','tem_ign_Seg','begS','endS','begS_half','endS_half','subject_id','count','lagvec','Numstages','NameStages','side_Tale','nE_id','userID','Expe','channels_to_reject')%% plot for subjects
end

%%
figure;
set(gcf,'Position',[440    70   509   728])
subplot(2,1,1);
for nSta=1:3
    gyt_PiratePlot(nSta-0.2,(mean_att(:,nSta)),0.2,0.5,'y',[0 0 1],'y',[],{'SizeData',36,'Marker','o'});
    gyt_PiratePlot(nSta+0.2,(mean_ign(:,nSta)),0.2,0.5,'y',[1 0 0],'y',[],{'SizeData',36,'Marker','o'});
end
set(gca,'XTick',1:3,'XTickLabel',{'Wake','NREM','REM'});
format_fig;
ylabel('Reconstruction score')

subplot(2,1,2);
for nSta=1:3
    gyt_PiratePlot(nSta,100*(mean_decoding(:,nSta)),0.4,0.5,'y',[1 1 1]*0.7,'y',[],{'SizeData',36,'Marker','o'});
end
set(gca,'YTick',0:20:100,'XTick',1:3,'XTickLabel',{'Wake','NREM','REM'})
ylim([0 110])
line(xlim,[1 1]*50,'Color','k','LineStyle','--')
format_fig;
ylabel('Decoding (%)')