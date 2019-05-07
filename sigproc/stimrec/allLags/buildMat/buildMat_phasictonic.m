NameStages={'tREM','pREM'};

phasic_subject(lenphases(2)>5*60)=[];%minimal 5 minutes
tonic_subject(lenphases(1)>5*60)=[];%minimal 5 minutes

rho_att_bySub_mean=zeros(length(subject_id),2);rho_ign_bySub_mean=zeros(length(subject_id),2);dec_bySub_mean=zeros(length(subject_id),2);
rho_att_bySub_Lag=zeros(length(subject_id),2,length(lagvec));rho_ign_bySub_Lag=zeros(length(subject_id),2,length(lagvec));dec_bySub_Lag_mean=zeros(length(subject_id),2,length(lagvec));
rho_att_bySub_Seg_mean=zeros(length(subject_id),2,length(begS));rho_ign_bySub_Seg_mean=zeros(length(subject_id),2,length(begS));dec_bySub_Seg_mean=zeros(length(subject_id),2,length(lagvec));
trialphasic=cell(1,length(subject_id));trialtonic=cell(1,length(subject_id));

for nS=1:length(subject_id)
    trialphasic{nS}=find(nrREMtrial{nS}>0);
    trialtonic{nS}=find(nrREMtrial{nS}==0);

    %% tonic
    nSta=1;
    rho_att_bySub_mean(nS,nSta)=mean(mean(tem_att{nS}(trialtonic{nS},:),2));
    rho_ign_bySub_mean(nS,nSta)=mean(mean(tem_ign{nS}(trialtonic{nS},:),2));
    dec_bySub_mean(nS,nSta)=100*mean(mean(tem_att{nS}(trialtonic{nS},:),2)>mean(tem_ign{nS}(trialtonic{nS},:),2));
    
    rho_att_bySub_Lag(nS,nSta,:)=mean(tem_att{nS}(trialtonic{nS},:),1);
    rho_ign_bySub_Lag(nS,nSta,:)=mean(tem_ign{nS}(trialtonic{nS},:),1);
    dec_bySub_Lag_mean(nS,nSta,:)=100*mean(tem_att{nS}(trialtonic{nS},:)>tem_ign{nS}(trialtonic{nS},:),1);
    for nSeg=1:length(begS)
        rho_att_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_att_Seg{nS}(trialtonic{nS},:,1),2)));
        rho_ign_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_ign_Seg{nS}(trialtonic{nS},:,2),2)));
        dec_bySub_Seg_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_Seg{nS}(trialtonic{nS},:,1),2)>mean(tem_ign{nS}(trialtonic{nS},:,1),2));
    end
    %% phasic
    nSta=2;
    rho_att_bySub_mean(nS,nSta)=mean(mean(tem_att{nS}(trialphasic{nS},:),2));
    rho_ign_bySub_mean(nS,nSta)=mean(mean(tem_ign{nS}(trialphasic{nS},:),2));
    dec_bySub_mean(nS,nSta)=100*mean(mean(tem_att{nS}(trialphasic{nS},:),2)>mean(tem_ign{nS}(trialphasic{nS},:),2));
    
    rho_att_bySub_Lag(nS,nSta,:)=mean(tem_att{nS}(trialphasic{nS},:),1);
    rho_ign_bySub_Lag(nS,nSta,:)=mean(tem_ign{nS}(trialphasic{nS},:),1);
    dec_bySub_Lag_mean(nS,nSta,:)=100*mean(tem_att{nS}(trialphasic{nS},:)>tem_ign{nS}(trialphasic{nS},:),1);
    for nSeg=1:length(begS)
        rho_att_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_att_Seg{nS}(trialphasic{nS},:,1),2)));
        rho_ign_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_ign_Seg{nS}(trialphasic{nS},:,2),2)));
        dec_bySub_Seg_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_Seg{nS}(trialtonic{nS},:,1),2)>mean(tem_ign_Seg{nS}(trialtonic{nS},:,1),2));
    end
    
end
goodsubj{1}=find(cellfun(@length,trialtonic)>=4);
goodsubj{2}=find(cellfun(@length,trialphasic)>=4);
