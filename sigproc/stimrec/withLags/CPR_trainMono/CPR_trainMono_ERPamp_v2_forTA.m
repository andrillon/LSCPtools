%% 

%%%%%%   !!! CPR_prepare_filtEEG must have been launched before and the files properly stored in the path

addpath('/home/lscp00/Programmes/matlab/toolbox/fieldtrip-20120521')
subject_id = {'003','004','005','006','008','009','010','011','012','013','014',...
    '015','016','017','019','020','021','022','023','024','026',...
    '027','028','029','030','031','032','033','034','035','036',...
    '037','038','039','040','041','042','043','044','045','047'};%
filterband=[2,8];%filtering
ExtractedDataPath=[];%where raw EEG and events are stored
DataPath=['?'];%where filtered data with eeg and audio are located
scoringPath=['?'];%where scoring files are stored
BehavPath=[];%where behavioral results files are stored

%% % parameters for ERP

NameStages={'Real','Jab','W','N2','REM'};
Numstages={[0,1],[0,1],[0,1],2,[5,6]};
maxdiff=[];
count=zeros(length(subject_id),length(Numstages));
att_meanERP=cell(length(subject_id),length(Numstages));
ign_meanERP=cell(length(subject_id),length(Numstages));
att_meanERP_half=cell(length(subject_id),length(Numstages));
ign_meanERP_half=cell(length(subject_id),length(Numstages));
side_Tale=cell(length(subject_id),length(Numstages));
nE_id=cell(length(subject_id),length(Numstages));
fs=500;
env_fs=100;%sample at 100Hz
% fs_resample=100;fltrate=1/fs_resample;subsamplerate=round(fs/fs_resample);
begS_half=[0,30];endS_half=[30,60];%time windows for 0-30s,30s-60s
channelSel=30;
filterERP=[0.5,30];
window_ERP=[-0.2,0.5];
baselineERP=[-0.2,0];
durminsilence=0.2;
% try
%     load([ExpePath,filesep,Expe 'artifacts']);
%     channels_to_reject=define_channel_artifacts(artifact_var);%
% catch
%     channels_to_reject=cell(1,length(subject_id));
% end
%%
for nS=1:length(subject_id)
    
    %%% Subj id
    subName=[Expe subject_id{nS}];
    SubID=subject_id{nS};
    fprintf('... Subject: %s\n',SubID);
    
    %%% Retrieve behavioral data
    SubjectBehavData=load([BehavPath filesep subName filesep 'Result_' subName '.mat']);
    fprintf('... ... behavioral data retrived\n')
    
    %%% Retrieve events EEG data
    load([ExtractedDataPath,'data_extracted' SubID '.mat'],'myevents','data_full')
    rerefData=data_full(channelSel,:)-mean(data_full([29,47],:),1);
    if filterERP(1)>0
        filt1=ft_preproc_highpassfilter(rerefData,500,filterERP(1),[],'firws','onepass-zerophase');
    else
        filt1=rerefData;
    end
    if filterERP(2)>0
        filt2=ft_preproc_lowpassfilter(rerefData,500,filterERP(2),[],'firws','onepass-zerophase');
    else
        filt2=filt1;
    end
                %%
    %%% retrieve structure of data : ordered by type of events, with
    %%% scoring in tagSleep
    CPR_trialevents
    CPR_trialstructure
    trial_sleep_scoring
    
    for st=1:4
        trialphase=trials_st(st);
        for nT=1:length(trialphase.id)
            nE=trialphase.id(nT);
            if st==1
                nSta=1;
            elseif st==2
                nSta=2;
            elseif st==3
                nSta=3;
            else
                nSta=find(cellfun(@(x) (trialphase.tagSleepScoring(nT)>=x(1) && trialphase.tagSleepScoring(nT)<=x(end)),Numstages(4:end)))+3;
                if isempty(nSta)
                    continue
                else
                    nSta=nSta(1);
                end
            end
            count(nS,nSta)=count(nS,nSta)+1;
            nE_id{nS,nSta}(count(nS,nSta))=nE;
            if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                side_Tale{nS,nSta}(count(nS,nSta))=1;
            elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                side_Tale{nS,nSta}(count(nS,nSta))=2;
            end
            
            saveName=sprintf('%s_T%02.0f_One_Pass_%g_%gHz_filtData',SubID,nE,filterband(1),filterband(2));
            load([DataPath filesep saveName '.mat'])
% filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
            
            % detect silences
            env1=filtData.env(:,1);
            true_timepoints = findSilences(env1,env_fs,'FIXED',0.025,durminsilence);
            silenceT1=true_timepoints(2,:);
            silenceT1(silenceT1/env_fs+window_ERP(1)<=0 | silenceT1/env_fs+window_ERP(2)>=length(env1)/env_fs)=[];%in seconds
            AvSilenceNumber{nS,nSta}(nE,1)=length(silenceT1);
            silenceT1=silenceT1/env_fs;
            
            erp1=[];torem=[];env_1=[];
            for nTsil=1:length(silenceT1)
                temp=filt2(round(trialphase.trials(nT,1)+silenceT1(nTsil)*fs+[window_ERP(1)*fs:window_ERP(2)*fs]))...
                    -mean(filt2(round(trialphase.trials(nT,1)+silenceT1(nTsil)*fs+[baselineERP(1)*fs:baselineERP(2)*fs])),2);
                %                 if sum(abs(temp)>erpartefact_limitamp)
%                     torem=[torem,nTsil];
%                     continue
%                 end
                erp1=[erp1 ; temp(1:end)];
%                 env=env1(silenceT1(nTsil)*env_fs+[window_ERP(1)*env_fs:window_ERP(2)*env_fs]);
%                 env_1=[env_1;env]
            end
            silenceT1(torem)=[];
            
            env2=filtData.env(:,2);
            true_timepoints = findSilences(env2,env_fs,'FIXED',0.025,durminsilence);
            silenceT2=true_timepoints(2,:);
            silenceT2(silenceT2/env_fs+window_ERP(1)<=0 | silenceT2/env_fs+window_ERP(2)>=length(env1)/env_fs)=[];%in seconds
            silenceT2=silenceT2/env_fs;
            AvSilenceNumber{nS,nSta}(nE,2)=length(silenceT2);
            
            if st>2
                erp2=[];torem=[];
                for nTsil=1:length(silenceT2)
                    temp=filt2(round(trialphase.trials(nT,1)+silenceT2(nTsil)*fs+[window_ERP(1)*fs:window_ERP(2)*fs]))...
                    -mean(filt2(round(trialphase.trials(nT,1)+silenceT2(nTsil)*fs+[baselineERP(1)*fs:baselineERP(2)*fs])),2);
%                     if sum(abs(temp)>100)
%                         torem=[torem,nTsil];
%                         continue
%                     end
                    erp2=[erp2 ; temp];
                end
                silenceT2(torem)=[];
            end
            
            if st>2
                if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                    att_meanERP{nS,nSta}(count(nS,nSta),:)= mean(erp1);
                    ign_meanERP{nS,nSta}(count(nS,nSta),:)= mean(erp2);
                elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                    att_meanERP{nS,nSta}(count(nS,nSta),:)= mean(erp2);
                    ign_meanERP{nS,nSta}(count(nS,nSta),:)= mean(erp1);
                end
            else
                att_meanERP{nS,nSta}(count(nS,nSta),:)= mean(erp1);
                ign_meanERP{nS,nSta}(count(nS,nSta),:)= mean(erp1);
            end
            
            %half
            if ~isempty(begS_half)
                for nSeg=1:length(begS_half)
                    silenceT1_Seg=find(silenceT1>begS_half(nSeg) & silenceT1<endS_half(nSeg));
                    silenceT2_Seg=find(silenceT2>begS_half(nSeg) & silenceT2<endS_half(nSeg));
                    AvSilenceNumber_Seg{nS,nSta}(count(nS,nSta),1,nSeg)=length(silenceT1_Seg);
                    AvSilenceNumber_Seg{nS,nSta}(count(nS,nSta),2,nSeg)=length(silenceT2_Seg);
                    if st>2
                        erp1_seg=erp1(silenceT1_Seg,:);
                        erp2_seg=erp2(silenceT2_Seg,:);
                        
                        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                            att_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(erp1_seg);
                            ign_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(erp2_seg);
                        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                            att_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(erp2_seg);
                            ign_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(erp1_seg);
                        end
                    else
                        erp1_seg=erp1(silenceT1_Seg,:);
                        att_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(erp1_seg);
                        ign_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(erp1_seg);
                        
                    end
                    
                end
                
            end
        end
    end
end

%% prepare plotting
goodsubj=cell(1,length(NameStages));
Rmatrix_diff_att_ign=[];Rmatrix_diff_att_ign_Seg=[];

att_meanERP_bySub=[];ign_meanERP_bySub=[];
att_meanERP_bySub_half=[];ign_meanERP_bySub_half=[];
ttestplot=1;

for nSta=1:length(NameStages)
    for nS=1:length(subject_id)
        if count(nS,nSta)>=4
            att_meanERP_bySub(nS,nSta,:)=nanmean(att_meanERP{nS,nSta},1);
            ign_meanERP_bySub(nS,nSta,:)=nanmean(ign_meanERP{nS,nSta},1);
                        
            if ~isempty(begS_half)
                for nSeg=1:length(begS_half)
                    att_meanERP_bySub_half(nS,nSta,nSeg,:)=nanmean(att_meanERP_half{nS,nSta}(:,:,nSeg),1);
                    ign_meanERP_bySub_half(nS,nSta,nSeg,:)=nanmean(ign_meanERP_half{nS,nSta}(:,:,nSeg),1);
                end
            end
            
        end
    end
    goodsubj{nSta}=find(count(:,nSta)>4);
end

%% plotting
figure; set(gcf,'position',[183        1103        1158         412])
timeP=[window_ERP(1):fltrate:window_ERP(2)];
color_Sta={[0,1,0],[0,0,1],[1,0,0],[1,1,0],[1,0,1],[0,1,1];[0,0.5,0],[0,0,0.5],[0.5,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]};
linew={'-','--'};
%%% reconstruction scores
for nSta=1:length(NameStages)
    subplot(1+~isempty(begS_half),length(NameStages),nSta)
    if ~isempty(goodsubj{nSta})
        hold on
        simpleTplot(timeP,{squeeze(att_meanERP_bySub(goodsubj{nSta},nSta,:)),squeeze(ign_meanERP_bySub(goodsubj{nSta},nSta,:))},0,cell2mat(color_Sta(1:2,nSta)),[0,0.05,0.05],linew{1},0.5,1,0);
        plot([timeP(1),timeP(end)],[0,0],'k--')
        plot([timeP(1),timeP(end)],[0,0],'k--')

        xlim(window_ERP)
        ylabel({'ERP on Cz','after silences in uV'})
        xlabel('time in s')
        set(gca,'FontSize',18,'FontWeight','bold')
        title([NameStages{nSta},', n=',num2str(length(goodsubj{nSta}))])
    end
end
%                     att_meanERP_bySub_half(nS,nSta,nSeg)=mean(att_meanERP_half{nS,nSta}(:,:,nSeg),1);

%%% half

for nSta=1:length(NameStages)
    subplot(1+~isempty(begS_half),length(NameStages),length(NameStages)+nSta)
    if ~isempty(goodsubj{nSta})
        hold on
        for nSeg=1:length(begS_half)
            simpleTplot(timeP,{squeeze(att_meanERP_bySub_half(goodsubj{nSta},nSta,nSeg,:)),squeeze(ign_meanERP_bySub_half(goodsubj{nSta},nSta,nSeg,:))},0,cell2mat(color_Sta(1:2,nSta)),[0,0.05,0.05],linew{nSeg},0.5,1,0);
        end
        plot([timeP(1),timeP(end)],[0,0],'k--')
        ylabel('ERP on Cz in uV')
        xlabel('time in s')
        set(gca,'FontSize',18,'FontWeight','bold')
        title({[NameStages{nSta},', n=',num2str(length(goodsubj{nSta}))],['-=1st half, --=2nd half']})
        xlim(window_ERP)
    end
end

