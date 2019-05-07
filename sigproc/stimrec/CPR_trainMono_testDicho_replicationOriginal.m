clear all
userID='Thomas';
filterband=[2,8];%filtering
Expe='CPR';
dateref='18-Sep-2018';%datestr(datetime('today'));%'26-Apr-2018';%datestr(datetime('today'));%%select results from which date to plot
lagvec=-50:0;% lag for linear model

if strcmp(userID,'Celia')
    CPR_path_config_Celia
elseif strcmp(userID,'Matthieu')
    pathGit=genpath('/home/lscp00/Work/EEGgit/Sommeil/CPR/');
    addpath(pathGit);
    CPR_path_define_Matthieu
elseif strcmp(userID,'Thomas')
    pathGit=('/Users/tand0009/WorkGit/projects/ext/Sommeil/CPR');
    addpath(genpath(pathGit));
    CPR_path_define_Thomas
end
NameStages={'Wake','N2','REM'};
subject_id = {'003','004','005','006','008','009','010','011','012','013','014',...
    '015','016','017','019','020','021','022','023','024','026',...
    '027','028','029','030','031','032','033','034','035','036',...
    '037','038','039','040','041','042','043','044','045','047'};%

BadChannels={'018',[23];
    '030',[23];
    '031',[23];
    '032',[53,54,55];
    '033',[23,24];
    '038',[23];
    '039',[23]};

%%
nS=0;
for nA=1:length(subject_id)
    
    %get subject id
    subName=[Expe subject_id{nA}];
    SubID=subject_id{nA};
    fprintf('... ... behavioural data retrieved\n')
    
    fprintf('... Subject: %s\n',SubID);
    SubjectBehavData=load([BehavPath filesep subName filesep 'Result_' subName '.mat']);
    
    %%% Retrieve events and EEG data
    saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,1,filterband(1),filterband(2));
    if exist([rawEEGPath,'data_extracted' SubID '.mat'])~=0 && exist([DataPath filesep saveName '.mat'])~=0
        load([rawEEGPath,'data_extracted' SubID '.mat'],'myevents')
    else
        continue;
    end
    nS=nS+1;
    %%% retrieve structure of data : ordered by type of events, with
    %%% scoring in tagSleep
    CPR_trialevents
    CPR_trialstructure
    %     trial_sleep_scoring
    
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
        load([DataPath filesep saveName '.mat']);
        
        train_data=[train_data ; filtData.eeg];
        train_stim=[train_stim ; filtData.env];
    end
    [g,rstim] = StimuliReconstruction (train_stim(:,1)', train_data, [], [], lagvec);
    
    fprintf('... ... test model on di-chotic\n')
    %% reconstruction Wake
    %Forced Wake Trials
    nSta=1;
    count(nS,nSta)=0;
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
        load([DataPath filesep saveName '.mat']);
        
        fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
        
        %stim reconstruction
        [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g, lagvec);
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
        mean_att(nS,nSta)=mean(tem_att{nS,nSta});
        mean_ign(nS,nSta)=mean(tem_ign{nS,nSta});
        mean_decoding(nS,nSta)=mean(tem_att{nS,nSta}>tem_ign{nS,nSta});
    end
    
    %Sleep Trials
    nSta=2;
    count(nS,nSta)=0;
    for nE=[trials_st(4).id]
        count(nS,nSta)=count(nS,nSta)+1;
        nE_id{nS,nSta}(count(nS,nSta))=nE;
        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
            side_Tale{nS,nSta}(count(nS,nSta))=1;
        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
            side_Tale{nS,nSta}(count(nS,nSta))=2;
        end
        %saveName for data
        saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        load([DataPath filesep saveName '.mat']);
        
        fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
        
        %stim reconstruction
        [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g, lagvec);
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
        mean_att(nS,nSta)=mean(tem_att{nS,nSta});
        mean_ign(nS,nSta)=mean(tem_ign{nS,nSta});
        mean_decoding(nS,nSta)=mean(tem_att{nS,nSta}>tem_ign{nS,nSta});
    end
    
    %Wake Sleep
    nSta=3;
    count(nS,nSta)=0;
    for nE=[trials_st(5).id]
        count(nS,nSta)=count(nS,nSta)+1;
        nE_id{nS,nSta}(count(nS,nSta))=nE;
        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
            side_Tale{nS,nSta}(count(nS,nSta))=1;
        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
            side_Tale{nS,nSta}(count(nS,nSta))=2;
        end
        %saveName for data
        saveName=sprintf('%s_T%02.0f_ButFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
        load([DataPath filesep saveName '.mat']);
        
        fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
        
        %stim reconstruction
        [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g, lagvec);
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
        mean_att(nS,nSta)=mean(tem_att{nS,nSta});
        mean_ign(nS,nSta)=mean(tem_ign{nS,nSta});
        mean_decoding(nS,nSta)=mean(tem_att{nS,nSta}>tem_ign{nS,nSta});
    end
    
%     %% training the model
%     clear g
%     fprintf('... ... training model on di-otic\n')
%     train_data=[];train_stim=[];
%     trials_training = [trials_st(1:2).id];
%     
%     for nE=trials_training
%         %saveName for data
%         saveName=sprintf('%s_T%02.0f_FIRFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
%         %load and reject artifacts if specified in channels_to_reject which
%         %channels to reject
%         load([DataPath filesep saveName '.mat']);
%         
%         train_data=[train_data ; filtData.eeg];
%         train_stim=[train_stim ; filtData.env];
%     end
%     [g,rstim] = StimuliReconstruction (train_stim(:,1)', train_data, [], [], lagvec);
%     
%     fprintf('... ... test model on di-chotic\n')
%     %% reconstruction Wake
%     %Forced Wake Trials
%     nSta=1;
%     count2(nS,nSta)=0;
%     for nE=[trials_st(3).id]
%         count2(nS,nSta)=count2(nS,nSta)+1;
%         nE_id2{nS,nSta}(count2(nS,nSta))=nE;
%         if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
%             side_Tale2{nS,nSta}(count2(nS,nSta))=1;
%         elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
%             side_Tale2{nS,nSta}(count2(nS,nSta))=2;
%         end
%         %saveName for data
%         saveName=sprintf('%s_T%02.0f_FIRFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
%         load([DataPath filesep saveName '.mat']);
%         
%         fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
%         
%         %stim reconstruction
%         [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g, lagvec);
%         [test_rho1, ~]=corr(rstimdum',filtData.env(:,1));
%         [test_rho2, ~]=corr(rstimdum',filtData.env(:,2));
%         
%         %take into account the side of the real speech
%         if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
%             tem_att2{nS,nSta}(count2(nS,nSta))=test_rho1;
%             tem_ign2{nS,nSta}(count2(nS,nSta))=test_rho2;
%         elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
%             tem_att2{nS,nSta}(count2(nS,nSta))=test_rho2;
%             tem_ign2{nS,nSta}(count2(nS,nSta))=test_rho1;
%         end
%         mean_att2(nS,nSta)=mean(tem_att2{nS,nSta});
%         mean_ign2(nS,nSta)=mean(tem_ign2{nS,nSta});
%         mean_decoding2(nS,nSta)=mean(tem_att2{nS,nSta}>tem_ign2{nS,nSta});
%     end
%     
%     %Sleep Trials
%     nSta=2;
%     count2(nS,nSta)=0;
%     for nE=[trials_st(4).id]
%         count2(nS,nSta)=count2(nS,nSta)+1;
%         nE_id2{nS,nSta}(count2(nS,nSta))=nE;
%         if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
%             side_Tale2{nS,nSta}(count2(nS,nSta))=1;
%         elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
%             side_Tale2{nS,nSta}(count2(nS,nSta))=2;
%         end
%         %saveName for data
%         saveName=sprintf('%s_T%02.0f_FIRFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
%         load([DataPath filesep saveName '.mat']);
%         
%         fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
%         
%         %stim reconstruction
%         [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g, lagvec);
%         [test_rho1, ~]=corr(rstimdum',filtData.env(:,1));
%         [test_rho2, ~]=corr(rstimdum',filtData.env(:,2));
%         
%         %take into account the side of the real speech
%         if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
%             tem_att2{nS,nSta}(count2(nS,nSta))=test_rho1;
%             tem_ign2{nS,nSta}(count2(nS,nSta))=test_rho2;
%         elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
%             tem_att2{nS,nSta}(count2(nS,nSta))=test_rho2;
%             tem_ign2{nS,nSta}(count2(nS,nSta))=test_rho1;
%         end
%         mean_att2(nS,nSta)=mean(tem_att2{nS,nSta});
%         mean_ign2(nS,nSta)=mean(tem_ign2{nS,nSta});
%         mean_decoding2(nS,nSta)=mean(tem_att2{nS,nSta}>tem_ign2{nS,nSta});
%     end
%     
%     %Wake Sleep
%     nSta=3;
%     count2(nS,nSta)=0;
%     for nE=[trials_st(5).id]
%         count2(nS,nSta)=count2(nS,nSta)+1;
%         nE_id2{nS,nSta}(count2(nS,nSta))=nE;
%         if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
%             side_Tale2{nS,nSta}(count2(nS,nSta))=1;
%         elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
%             side_Tale2{nS,nSta}(count2(nS,nSta))=2;
%         end
%         %saveName for data
%         saveName=sprintf('%s_T%02.0f_FIRFilt_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
%         load([DataPath filesep saveName '.mat']);
%         
%         fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
%         
%         %stim reconstruction
%         [gdum,rstimdum] = StimuliReconstruction ([], [], filtData.eeg, g, lagvec);
%         [test_rho1, ~]=corr(rstimdum',filtData.env(:,1));
%         [test_rho2, ~]=corr(rstimdum',filtData.env(:,2));
%         
%         %take into account the side of the real speech
%         if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
%             tem_att2{nS,nSta}(count2(nS,nSta))=test_rho1;
%             tem_ign2{nS,nSta}(count2(nS,nSta))=test_rho2;
%         elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
%             tem_att2{nS,nSta}(count2(nS,nSta))=test_rho2;
%             tem_ign2{nS,nSta}(count2(nS,nSta))=test_rho1;
%         end
%         mean_att2(nS,nSta)=mean(tem_att2{nS,nSta});
%         mean_ign2(nS,nSta)=mean(tem_ign2{nS,nSta});
%         mean_decoding2(nS,nSta)=mean(tem_att2{nS,nSta}>tem_ign2{nS,nSta});
%     end
end

%%
figure;
set(gcf,'Position',[440    70   509   728])
subplot(2,1,1);
for nSta=1:3
    gyt_PiratePlot(nSta-0.2,(mean_att(:,nSta)),0.2,0.5,'y',[0 0 1],'y',[],{'SizeData',36,'Marker','o'});
    gyt_PiratePlot(nSta+0.2,(mean_ign(:,nSta)),0.2,0.5,'y',[1 0 0],'y',[],{'SizeData',36,'Marker','o'});
end
set(gca,'XTick',1:3,'XTickLabel',{'Forced Wake','Sleep Trials','Wake Sleep'});
format_fig;
ylabel('Reconstruction score')

subplot(2,1,2);
for nSta=1:3
    gyt_PiratePlot(nSta,100*(mean_decoding(:,nSta)),0.4,0.5,'y',[1 1 1]*0.7,'y',[],{'SizeData',36,'Marker','o'});
end
set(gca,'YTick',0:20:100,'XTick',1:3,'XTickLabel',{'Forced W','Sleep Tr','Wake Sleep'})
ylim([0 110])
line(xlim,[1 1]*50,'Color','k','LineStyle','--')
format_fig;
ylabel('Decoding (%)')


% %%
% figure;
% set(gcf,'Position',[440    70   509   728])
% subplot(2,1,1);
% for nSta=1:3
%     gyt_PiratePlot(nSta-0.2,(mean_att2(:,nSta)),0.2,0.5,'y',[0 0 1],'y',[],{'SizeData',36,'Marker','o'});
%     gyt_PiratePlot(nSta+0.2,(mean_ign2(:,nSta)),0.2,0.5,'y',[1 0 0],'y',[],{'SizeData',36,'Marker','o'});
% end
% set(gca,'XTick',1:3,'XTickLabel',{'Forced W','Sleep Tr','Wake Sleep'});
% format_fig;
% ylabel('Reconstruction score')
% 
% subplot(2,1,2);
% for nSta=1:3
%     gyt_PiratePlot(nSta,100*(mean_decoding2(:,nSta)),0.4,0.5,'y',[1 1 1]*0.7,'y',[],{'SizeData',36,'Marker','o'});
% end
% set(gca,'YTick',0:20:100,'XTick',1:3,'XTickLabel',{'Forced W','Sleep Tr','Wake Sleep'})
% ylim([0 110])
% line(xlim,[1 1]*50,'Color','k','LineStyle','--')
% format_fig;
% ylabel('Decoding (%)')