NameStages={'tREM','pREM'};
rho_att_bySub_mean=zeros(length(subject_id),2);rho_ign_bySub_mean=zeros(length(subject_id),2);dec_bySub_mean=zeros(length(subject_id),2);
rho_att_bySub_Seg_mean=zeros(length(subject_id),2,length(begS));rho_ign_bySub_Seg_mean=zeros(length(subject_id),2,length(begS));dec_bySub_Seg_mean=zeros(length(subject_id),2,length(lagvec));
rho_att_bySub_half_mean=zeros(length(subject_id),2,length(begS));rho_ign_bySub_half_mean=zeros(length(subject_id),2,length(begS));dec_bySub_half_mean=zeros(length(subject_id),2,length(lagvec));
trialphasic=cell(1,length(subject_id));trialtonic=cell(1,length(subject_id));
retrictLags=[0:50];restrictLagsvec=find(ismember(lagvec,retrictLags));
lagsplot=lagvec(restrictLagsvec);
rho_att_bySub_Lag=zeros(length(subject_id),2,length(restrictLagsvec));rho_ign_bySub_Lag=zeros(length(subject_id),2,length(restrictLagsvec));dec_bySub_Lag_mean=zeros(length(subject_id),2,length(restrictLagsvec));
goodsubj=cell(1,2);

for nS=1:length(subject_id)
    trialphasic{nS}=find(nrREMtrial{nS}>0);
    trialtonic{nS}=find(nrREMtrial{nS}==0);
    if ~isempty(trialtonic{nS})
    %% tonic
    nSta=1;
    rho_att_bySub_mean(nS,nSta)=mean(mean(tem_att{nS}(trialtonic{nS},restrictLagsvec),2));
    rho_ign_bySub_mean(nS,nSta)=mean(mean(tem_ign{nS}(trialtonic{nS},restrictLagsvec),2));
    dec_bySub_mean(nS,nSta)=100*mean(mean(tem_att{nS}(trialtonic{nS},restrictLagsvec),2)>mean(tem_ign{nS}(trialtonic{nS},restrictLagsvec),2));
    
    rho_att_bySub_Lag(nS,nSta,:)=mean(tem_att{nS}(trialtonic{nS},restrictLagsvec),1);
    rho_ign_bySub_Lag(nS,nSta,:)=mean(tem_ign{nS}(trialtonic{nS},restrictLagsvec),1);
    dec_bySub_Lag_mean(nS,nSta,:)=100*mean(tem_att{nS}(trialtonic{nS},restrictLagsvec)>tem_ign{nS}(trialtonic{nS},restrictLagsvec),1);
    if ~isempty(begS_half)
        for nSeg=1:length(begS_half)
            rho_att_bySub_half_mean(nS,nSta,nSeg)=mean(mean(tem_att_half{nS}(:,restrictLagsvec,nSeg),2));
            rho_ign_bySub_half_mean(nS,nSta,nSeg)=mean(mean(tem_ign_half{nS}(:,restrictLagsvec,nSeg),2));
            dec_bySub_half_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_half{nS}(:,restrictLagsvec,nSeg),2)>mean(tem_ign_half{nS}(:,restrictLagsvec,nSeg),2));
        end
    end
    if ~isempty(begS)
        for nSeg=1:length(begS)
            rho_att_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_att_Seg{nS}(trialtonic{nS},restrictLagsvec,nSeg),2)));
            rho_ign_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_ign_Seg{nS}(trialtonic{nS},restrictLagsvec,nSeg),2)));
            dec_bySub_Seg_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_Seg{nS}(trialtonic{nS},restrictLagsvec,nSeg),2)>mean(tem_ign_Seg{nS}(trialtonic{nS},restrictLagsvec,nSeg),2));
        end
    end
    end
    %% phasic
    if ~isempty(trialphasic{nS})
    nSta=2;
    rho_att_bySub_mean(nS,nSta)=mean(mean(tem_att{nS}(trialphasic{nS},restrictLagsvec),2));
    rho_ign_bySub_mean(nS,nSta)=mean(mean(tem_ign{nS}(trialphasic{nS},restrictLagsvec),2));
    dec_bySub_mean(nS,nSta)=100*mean(mean(tem_att{nS}(trialphasic{nS},restrictLagsvec),2)>mean(tem_ign{nS}(trialphasic{nS},restrictLagsvec),2));
    
    rho_att_bySub_Lag(nS,nSta,:)=mean(tem_att{nS}(trialphasic{nS},restrictLagsvec),1);
    rho_ign_bySub_Lag(nS,nSta,:)=mean(tem_ign{nS}(trialphasic{nS},restrictLagsvec),1);
    dec_bySub_Lag_mean(nS,nSta,:)=100*mean(tem_att{nS}(trialphasic{nS},restrictLagsvec)>tem_ign{nS}(trialphasic{nS},restrictLagsvec),1);
    if ~isempty(begS_half)
        for nSeg=1:length(begS_half)
            rho_att_bySub_half_mean(nS,nSta,nSeg)=mean(mean(tem_att_half{nS}(:,restrictLagsvec,nSeg),2));
            rho_ign_bySub_half_mean(nS,nSta,nSeg)=mean(mean(tem_ign_half{nS}(:,restrictLagsvec,nSeg),2));
            dec_bySub_half_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_half{nS}(:,restrictLagsvec,nSeg),2)>mean(tem_ign_half{nS}(:,restrictLagsvec,nSeg),2));
        end
    end
    if ~isempty(begS)
        for nSeg=1:length(begS)
            rho_att_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_att_Seg{nS}(trialphasic{nS},restrictLagsvec,nSeg),2)));
            rho_ign_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_ign_Seg{nS}(trialphasic{nS},restrictLagsvec,nSeg),2)));
            dec_bySub_Seg_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_Seg{nS}(trialphasic{nS},restrictLagsvec,nSeg),2)>mean(tem_ign_Seg{nS}(trialphasic{nS},restrictLagsvec,nSeg),2));
        end
    end
    end
end
goodsubj{1}=find(cellfun(@length,trialtonic)>=4);
goodsubj{2}=find(cellfun(@length,trialphasic)>=4);
