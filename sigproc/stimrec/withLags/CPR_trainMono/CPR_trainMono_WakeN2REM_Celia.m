%% initialisation
subject_id = {'006','008','009','010','011','012','013','014',...
    '015','016','017','019','020','021','022','023','024','026',...
    '027','028','029','030','031','032','033','034','035','036',...
    '037','038','039'};%'011',;%,'050'};%,'090','100','130','150','160','170','180','200','210','220'};%,'190'};

Expe = 'CPR';
filterband=[0.5,8];%filtering
userID='Celia';

%%% path define
if strcmp(userID,'Celia')
    CPR_path_config_Celia
elseif strcmp(userID,'Matthieu')
    CPR_path_define_Matthieu
elseif strcmp(userID,'Thomas')
    CPR_path_define_Thomas
end


%%% parameters for reconstruction
lagvec=[0:50];% lag for linear model
fltrate=0.01;%sampled at 10Hz
begS=[1,30./fltrate];endS=[30/fltrate+1,60/fltrate];%time windows for 0-30s,30s-60s

NameStages={'W', 'N2','REM'};
Numstages={[0,1],2,[5,6]};
maxdiff=[];
count=zeros(length(subject_id),length(Numstages));
tem_att=cell(length(subject_id),length(Numstages));
tem_ign=cell(length(subject_id),length(Numstages));

%% Reconstruction
for nS=1:length(subject_id)
    
    subName=[Expe subject_id{nS}];
    SubID=subName(4:end);
    fprintf('... Subject: %s\n',SubID);
    
    %%% Retrieve behavioural data
    SubjectBehavData=load([BehavPath filesep subName filesep 'Result_' subName '.mat']);
    load([findPath filesep 'trials_str' SubID '.mat']);
    fprintf('... ... behavioral data retrieved\n')
    
    %%% Retrieve events EEG data
        load([rawEEGPath,'data_extracted' SubID '.mat'],'myevents')
        hdr=ft_read_event([rawPath,filesep subName,'.raw']);
    
    
    %%% retrieve structure of data : ordered by type of events, with
    %%% scoring in tagSleep
    CPR_trialevents
    CPR_trialstructure
    trial_sleep_scoring
    
    clear g
    %%% training the model
    fprintf('... ... training model on di-otic\n')
    train_data=[];train_stim=[];
    trials_training = [trials_st(1:2).id];
    for nE=trials_training

        if strcmp(userID,'Celia');saveName=sprintf('%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        elseif strcmp(userID,'Matthieu');saveName=sprintf('%s_T%02.0f_filtData',SubID,nE);end
        load([DataPath filesep saveName '.mat']);        
        train_data=[train_data ; filtData.eeg];
        train_stim=[train_stim ; filtData.env];
    end
    for nLag=1:length(lagvec)
        [g(:,nLag),rstim] = StimuliReconstruction (train_stim(:,1)', train_data, [], [], -lagvec(nLag));
    end
    
    fprintf('... ... test model on di-chotic\n')
    nSta=1;
    
    %Forced Wake Trials
    for nE=[trials_st(3).id]
        %       fprintf('.. %g/40 ..\n',nE)
        
        if strcmp(userID,'Celia');saveName=sprintf('%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        elseif strcmp(userID,'Matthieu');saveName=sprintf('%s_T%02.0f_filtData',SubID,nE);end
        load([DataPath filesep saveName '.mat']);
        %         saveName=sprintf('%s_T%02.0f_filtData',SubID,nE);
        %         load([filteredEEGPath filesep saveName])
        
        count(nS,nSta)=count(nS,nSta)+1;
        nE_id{nS,nSta}(count(nS,nSta))=nE;
        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
            side_Tale_Sleep{nS,nSta}(count(nS,nSta))=1;
        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
            side_Tale_Sleep{nS,nSta}(count(nS,nSta))=2;
        end
        
        fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
        for nLag=1:length(lagvec)
            
            [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g(:,nLag), -lagvec(nLag));
            [test_rho1, ~]=corr(rstimdum',filtData.env(:,1));
            [test_rho2, ~]=corr(rstimdum',filtData.env(:,2));
            
            if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                tem_att{nS,nSta}(count(nS,nSta),nLag)=test_rho1;
                tem_ign{nS,nSta}(count(nS,nSta),nLag)=test_rho2;
            elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                tem_att{nS,nSta}(count(nS,nSta),nLag)=test_rho2;
                tem_ign{nS,nSta}(count(nS,nSta),nLag)=test_rho1;
            end
        end
        
            %30s-60s
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
    
    
    %Sleeptest trials
    for nT=1:length([trials_st(4).id])
        
        %select N2 and REM
        nSta=find(cellfun(@(x) ismember(trials_st(4).tagSleepScoring(nT),x),Numstages(2:end)))+1;
        if isempty(nSta)
            continue
        end
        
        if strcmp(userID,'Celia');saveName=sprintf('%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        elseif strcmp(userID,'Matthieu');saveName=sprintf('%s_T%02.0f_filtData',SubID,nE);end
        load([DataPath filesep saveName '.mat']);
        
        count(nS,nSta)=count(nS,nSta)+1;
        nE=trials_st(4).id(nT);
        nE_id{nS,nSta}(count(nS,nSta))=nE;
        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
            side_Tale_Sleep{nS,nSta}(count(nS,nSta))=1;
        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
            side_Tale_Sleep{nS,nSta}(count(nS,nSta))=2;
        end
        
        fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
        for nLag=1:length(lagvec)
            [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g(:,nLag), -lagvec(nLag));
            
            [test_rho1, ~]=corr(rstimdum',filtData.env(:,1));
            [test_rho2, ~]=corr(rstimdum',filtData.env(:,2));
            
            %             %mTRF
            %             [mTRF1,t] = mTRFtrain(filtData.env(thisSeg,1),filtData.eeg(thisSeg,:),100,0,-150,450,100);%stim,EEG,fs,map,tmin,tmax,lambda
            %             [mTRF2,t] = mTRFtrain(filtData.env(thisSeg,2),filtData.eeg(thisSeg,:),100,0,-150,450,100);%stim,EEG,fs,map,tmin,tmax,lambda
            
            if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                tem_att{nS,nSta}(count(nS,nSta),nLag)=test_rho1;
                tem_ign{nS,nSta}(count(nS,nSta),nLag)=test_rho2;
                %                 TRF_att{nSta}(count(nSta),1:length(t))=mTRF1;
                %                 TRF_ign{nSta}(count(nSta),1:length(t))=mTRF2;
            elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                tem_att{nS,nSta}(count(nS,nSta),nLag)=test_rho2;
                tem_ign{nS,nSta}(count(nS,nSta),nLag)=test_rho1;
                %                 TRF_att{nSta}(count(nSta),1:length(t))=mTRF2;
                %                 TRF_ign{nSta}(count(nSta),1:length(t))=mTRF1;
            end
        end
        
        %30s-60s
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
end
a=datetime;
save([savePath,filesep,'Reconstruction',num2str(filterband(1)),'to',num2str(filterband(2)),datestr(a) '.mat'],'tem_att','tem_ign','tem_att_Seg','tem_ign_Seg','begS','endS','subject_id','count','lagvec','Numstages','NameStages','side_Tale_Sleep','nE_id')%% plot for subjects

%%% plot Results

%%

goodsubj=cell(1,length(Numstages));
for nSta=1:length(Numstages)
    for nS=1:length(subject_id)
        if count(nS,nSta)>=4
            rho_att_bySub_mean(nS,nSta)=mean(mean(tem_att{nS,nSta},2));
            rho_ign_bySub_mean(nS,nSta)=mean(mean(tem_ign{nS,nSta},2));
            dec_bySub_mean(nS,nSta,:)=100*mean(mean(tem_att{nS,nSta},2)>mean(tem_ign{nS,nSta},2));
            
            rho_att_bySub_Lag(nS,nSta,:)=mean(tem_att{nS,nSta},1);
            rho_ign_bySub_Lag(nS,nSta,:)=mean(tem_ign{nS,nSta},1);
            dec_bySub_Lag_mean(nS,nSta,:)=100*mean(tem_att{nS,nSta}>tem_ign{nS,nSta},1);
            
            if ~isempty(begS)
                for nSeg=1:length(begS)
                    rho_att_bySub_Seg_mean(nS,nSta,nSeg)=mean(mean(tem_att_Seg{nS,nSta}(:,:,nSeg),2));
                    rho_ign_bySub_Seg_mean(nS,nSta,nSeg)=mean(mean(tem_ign_Seg{nS,nSta}(:,:,nSeg),2));
                    dec_bySub_Seg_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_Seg{nS,nSta}(:,:,nSeg),2)>mean(tem_ign_Seg{nS,nSta}(:,:,nSeg),2));
                end
            end
        end
    end
    goodsubj{nSta}=find(count(:,nSta)>4);
end
%% plot figure1
plotFigure1_Celia

% %% plot trimadere
% plot_Figure1Lag