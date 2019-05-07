%% initialisation
subject_id = {'006','008','009','010','011','012','013','014',...
    '015','016','017','019','020','021','022','023','024','026',...
    '027','028','029','030','031','032','033','034','035','036','037','038','039'};%'011',;%,'050'};%,'090','100','130','150','160','170','180','200','210','220'};%,'190'};

Expe = 'CPR';
filterband=[0.5,8];%filtering
userID='Matthieu';

%%% path define
if strcmp(userID,'Celia')
    CPR_path_config_Celia
elseif strcmp(userID,'Matthieu')
    pathGit =genpath('/home/lscp00/Work/EEGgit/Sommeil/');
    addpath(pathGit);
    CPR_path_define_Matthieu
elseif strcmp(userID,'Thomas')
    CPR_path_define_Thomas
end

% matthieuPath='/media/lscp00/My Passport/';
% ExpePath=[matthieuPath Expe filesep]

lagvec=[0:50];% lag for linear model
fltrate=0.01;%sampled at 10Hz
% begS=[1,30./fltrate];endS=[30/fltrate+1,60/fltrate];%time windows for 0-30s,30s-60s
maxdiff=[];
NumStages={[0,1],[0,1]};
NameStages={'real','jab'};
tem_mono=cell(length(subject_id),2);
maxdiff=[];

count=zeros(length(subject_id),length(NameStages));
% tem_mono=cell(length(subject_id),length(NameStages));
load([ExpePath,filesep,Expe 'artifacts']);
channels_to_reject=cell(1,length(subject_id));%define_channel_artifacts(artifact_var);
% cell(1,length(subject_id));%
%% Reconstruction
% 
%     13
%     22
%     25
%     30
for nS=length(tem_mono1)%length(subject_id)%setdiff(19:length(subject_id),[13,22,25,30,15,16]);%
    
    %%% Import in SPM
    subName=[Expe subject_id{nS}];
    SubID=subject_id{nS};
    fprintf('... Subject: %s\n',SubID);
    
    %%% Retrive behavioral data
    SubjectBehavData=load([BehavPath filesep subName filesep 'Result_' subName '.mat']);
    fprintf('... ... behavioral data retrived\n')
    
    %%% Retrieve events EEG data
    load([ExtractedDataPath,'data_extracted' SubID '.mat'],'myevents');
    %%% retrieve structure of data : ordered by type of events, with
    %%% scoring in tagSleep
%     CPR_trialevents
%     CPR_trialstructure
%     trial_sleep_scoring
    saveName=['trials_str',num2str(SubID)];
    load([trialstPath filesep saveName],'trials_st')
     
    %%% Epoching and saving
    train_data=[];train_stim=[];
    trials_training = [trials_st(1:2).id];
    %do the leave-one out
    %%
    for nr=7:length(trials_training)
        
        %select and arrange trial for leave-one out
        nE=trials_training(nr);
        trials_id=trials_training(setdiff(1:length(trials_training),nE));
        g=[];
        if ismember(nE,trials_st(1).id)
            attorign=1;
        elseif  ismember(nE,trials_st(2).id)
            attorign=2;
        end
        nrtype=find(trials_st(attorign).id==nE);
        
        %train on all other trials
        fprintf('... ... auto-test model\n')
        for nE=trials_id
            fprintf('.. %g/40 ..\n',nE)
            clear filtData
            
            if strcmp(userID,'celia');saveName=sprintf(['%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2)]);
            elseif strcmp(userID,'Matthieu');saveName=sprintf('%s_T%02.0f_filtData',SubID,nE);
            elseif strcmp(userID,'Thomas');saveName=sprintf('%s_T%02.0f_filtData',SubID,nE);end
            filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
            %             saveName=sprintf('%s_T%02.0f_filtData.mat',SubID,nE);
            %             load([DataPath filesep saveName]);
            filtData
            train_data=[train_data ; filtData.eeg];
            train_stim=[train_stim ; filtData.env];
        end
        for nLag=1:length(lagvec)
            [g(:,nLag),rstim] = StimuliReconstruction(train_stim(:,1)', train_data, [], [], -lagvec(nLag));
        end
        fprintf('... ... test model on diotic\n')
        
        %select the trial to test
        nE=trials_training(nr);
        if strcmp(userID,'celia');saveName=sprintf(['%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2)]);
        elseif strcmp(userID,'Matthieu');saveName=sprintf('%s_T%02.0f_filtData',SubID,nE);end
        filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
        %         saveName=sprintf('%s_T%02.0f_filtData',SubID,nE);
        %         load([DataPath filesep saveName])
        for nLag=1:length(lagvec)
            [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g(:,nLag), -nLag);
            [test_rho, ~]=corr(rstimdum',filtData.env(:,1));
            tem_mono1{nS,attorign}(nrtype,nLag)=test_rho;
        end
    end
    
end
%%
a=datetime;
save([savePath,filesep,userID,'_ControlLeaveoneOut',num2str(filterband(1)),'to',num2str(filterband(2)),'.mat' ],'tem_mono1','tem_mono','channels_to_reject','subject_id','count','lagvec','NumStages','NameStages')%% plot for subjects

%%% plot Results

%%
rho_mono_bySub_mean=zeros(length(subject_id),2);
for nS=1:length(subject_id)
    for type=1:2
        try
            rho_mono_bySub_mean(nS,type)=mean(mean(tem_mono{nS,type},2));
            rho_att_bySub_Lag(nS,type,1:51)=mean(tem_mono{nS,type},1);
        end
    end
end
%%
goodsubj=find(cellfun(@isempty,tem_mono(:,1))==0);
figure;
phase_ind=[1,2];colorplot=['b','r'];
hold on;
for type=1:2
    bar(phase_ind(type),mean(rho_mono_bySub_mean(goodsubj,type)),colorplot(type))
    errorbar(phase_ind(type),mean(rho_mono_bySub_mean(goodsubj,type)),std(rho_mono_bySub_mean(goodsubj,type))/sqrt(length(goodsubj)),colorplot(type))
    
    [h,p]=ttest_plot(gcf,gca,phase_ind(type),rho_mono_bySub_mean(goodsubj,type))
end
[h,p]=ttest_plot(gcf,gca,phase_ind,rho_mono_bySub_mean(goodsubj,1),rho_mono_bySub_mean(goodsubj,2))

NameStages={'real','jabb'};
plot([0.5,2.5],[0,0],'k--')
set(gca,'XTick',1:2,'XTickLabel',NameStages)
ylabel('reconstruction score')
set(gca,'FontSize',18,'FontWeight','bold')
title(['rec.scores, blue=real, red=jab,n=',length(goodsubj)])
