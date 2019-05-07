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
erpartefact_limitamp=100;
% try
%     load([ExpePath,filesep,Expe 'artifacts']);
%     channels_to_reject=define_channel_artifacts(artifact_var);%
% catch
    channels_to_reject=cell(1,length(subject_id));
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
    load([ExtractedDataPath,'data_extracted' SubID '.mat'],'myevents')
%     rerefData=data_full(channelSel,:)-mean(data_full([29,47],:),1);
%     if filterERP(1)>0
%         filt1=ft_preproc_highpassfilter(rerefData,500,filterERP(1),[],'firws','onepass-zerophase');
%     else
%         filt1=rerefData;
%     end
%     if filterERP(2)>0
%         filt2=ft_preproc_lowpassfilter(rerefData,500,filterERP(2),[],'firws','onepass-zerophase');
%     else
%         filt2=filt1;
%     end
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
            filtData=prepare_filtData([DataPath filesep saveName '.mat'],channels_to_reject{nS});
            
            % detect silences
            env1=filtData.env(:,1);
            silenceT=find(env1<0.025);
            sep_silenceT=diff(silenceT);
            silenceT1=silenceT([sep_silenceT ; 1]./env_fs>durminsilence);
            silenceT1(silenceT1/env_fs+window_ERP(1)<=0 | silenceT1/env_fs+window_ERP(2)>=length(env1)/env_fs)=[];%in seconds
            AvSilenceNumber{nS,nSta}(nE,1)=length(silenceT1);
            silenceT1=silenceT1/env_fs;
            
            env1_all=[];torem=[];env_1=[];
            for nTsil=1:length(silenceT1)
                %                 temp=filt2(round(trialphase.trials(nT,1)+silenceT1(nTsil)*fs+[window_ERP(1)*fs:window_ERP(2)*fs]))...
                %                     -mean(filt2(round(trialphase.trials(nT,1)+silenceT1(nTsil)*fs+[baselineERP(1)*fs:baselineERP(2)*fs])),2);
                %                 %                 if sum(abs(temp)>erpartefact_limitamp)
                % %                     torem=[torem,nTsil];
                % %                     continue
                % %                 end
                %                 erp1=[erp1 ; temp(1:end)];
                env=env1(round(silenceT1(nTsil)*env_fs+[window_ERP(1)*env_fs:window_ERP(2)*env_fs]))';
                env1_all=[env1_all;env];
            end
            silenceT1(torem)=[];
            
            env2=filtData.env(:,2);
            silenceT=find(env2<0.025);
            sep_silenceT=diff(silenceT);
            silenceT2=silenceT([sep_silenceT ; 1]./env_fs>durminsilence);
            silenceT2(silenceT2/env_fs+window_ERP(1)<=0 | silenceT2/env_fs+window_ERP(2)>=length(env1)/env_fs)=[];%in seconds
            silenceT2=silenceT2/env_fs;
            AvSilenceNumber{nS,nSta}(nE,2)=length(silenceT2);
            
            if st>2
                env2_all=[];torem=[];
                for nTsil=1:length(silenceT2)
                    %                     temp=filt2(round(trialphase.trials(nT,1)+silenceT2(nTsil)*fs+[window_ERP(1)*fs:window_ERP(2)*fs]))...
                    %                     -mean(filt2(round(trialphase.trials(nT,1)+silenceT2(nTsil)*fs+[baselineERP(1)*fs:baselineERP(2)*fs])),2);
                    %                     if sum(abs(temp)>100)
                    %                         torem=[torem,nTsil];
                    %                         continue
                    %                     end
                    env2_all=[env2_all ; env2(round(silenceT2(nTsil)*env_fs+[window_ERP(1)*env_fs:window_ERP(2)*env_fs]))'];
                end
                %                 silenceT2(torem)=[];
            end
            
            if st>2
                if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                    att_meanERP{nS,nSta}(count(nS,nSta),:)= mean(env1_all);
                    ign_meanERP{nS,nSta}(count(nS,nSta),:)= mean(env2_all);
                elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                    att_meanERP{nS,nSta}(count(nS,nSta),:)= mean(env2_all);
                    ign_meanERP{nS,nSta}(count(nS,nSta),:)= mean(env1_all);
                end
            else
                att_meanERP{nS,nSta}(count(nS,nSta),:)= mean(env1_all);
                ign_meanERP{nS,nSta}(count(nS,nSta),:)= mean(env1_all);
            end
            
            %half
            if ~isempty(begS_half)
                for nSeg=1:length(begS_half)
                    silenceT1_Seg=find(silenceT1>begS_half(nSeg) & silenceT1<endS_half(nSeg));
                    silenceT2_Seg=find(silenceT2>begS_half(nSeg) & silenceT2<endS_half(nSeg));
                    AvSilenceNumber_Seg{nS,nSta}(count(nS,nSta),1,nSeg)=length(silenceT1_Seg);
                    AvSilenceNumber_Seg{nS,nSta}(count(nS,nSta),2,nSeg)=length(silenceT2_Seg);
                    if st>2
                        env1_seg=env1_all(silenceT1_Seg,:);
                        env2_seg=env2_all(silenceT2_Seg,:);
                        
                        if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                            att_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(env1_seg);
                            ign_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(env2_seg);
                        elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                            att_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(env2_seg);
                            ign_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(env1_seg);
                        end
                    else
                        env1_seg=env1_all(silenceT1_Seg,:);
                        att_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(env1_seg);
                        ign_meanERP_half{nS,nSta}(count(nS,nSta),:,nSeg)=mean(env1_seg);
                        
                    end
                    
                end
                
            end
        end
    end
end