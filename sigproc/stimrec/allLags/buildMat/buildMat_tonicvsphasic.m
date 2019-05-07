NameStages={'tREM','pREM'};
clear goodsubj
goodsubj{1}=find(lenphases(:,1)>5*60);%minimal 5 minutes
goodsubj{2}=find(lenphases(:,2)>5*60);%minimal 5 minutes

rho_att_bySub_mean=zeros(length(subject_id),2);rho_ign_bySub_mean=zeros(length(subject_id),2);dec_bySub_mean=zeros(length(subject_id),2);
rho_att_bySub_Lag=zeros(length(subject_id),2,length(lagvec));rho_ign_bySub_Lag=zeros(length(subject_id),2,length(lagvec));dec_bySub_Lag_mean=zeros(length(subject_id),2,length(lagvec));
trialphasic=cell(1,length(subject_id));trialtonic=cell(1,length(subject_id));
begS=[];begS_half=[];
retrictLags=[0:50];restrictLagsvec=find(ismember(lagvec,retrictLags));
lagsplot=lagvec(restrictLagsvec);

for nS=1:length(subject_id)
    
    
    %% tonic
    nSta=1;
    rho_att_bySub_mean(nS,nSta)=mean(rho_real_tonic(nS,restrictLagsvec));
    rho_ign_bySub_mean(nS,nSta)=mean(rho_jab_tonic(nS,restrictLagsvec));
    dec_bySub_mean(nS,nSta)=100*(rho_att_bySub_mean(nS,nSta)>rho_ign_bySub_mean(nS,nSta));
    
    rho_att_bySub_Lag(nS,nSta,:)=rho_real_tonic(nS,restrictLagsvec);
    rho_ign_bySub_Lag(nS,nSta,:)=rho_jab_tonic(nS,restrictLagsvec);
    dec_bySub_Lag_mean(nS,nSta,:)=100*(rho_real_tonic(nS,restrictLagsvec)>rho_jab_tonic(nS,restrictLagsvec));
    %     for nSeg=1:length(begS)
    %         rho_att_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_att_Seg{nS}(trialtonic{nS},:,1),2)));
    %         rho_ign_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_ign_Seg{nS}(trialtonic{nS},:,2),2)));
    %         dec_bySub_Seg_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_Seg{nS}(trialtonic{nS},:,1),2)>mean(tem_ign{nS}(trialtonic{nS},:,1),2));
    %     end
    %% phasic
    nSta=2;
    rho_att_bySub_mean(nS,nSta)=mean(rho_real_phasic(nS,restrictLagsvec));
    rho_ign_bySub_mean(nS,nSta)=mean(rho_jab_phasic(nS,restrictLagsvec));
    dec_bySub_mean(nS,nSta)=100*(rho_att_bySub_mean(nS,nSta)>rho_ign_bySub_mean(nS,nSta));
    
    rho_att_bySub_Lag(nS,nSta,:)=rho_real_phasic(nS,restrictLagsvec);
    rho_ign_bySub_Lag(nS,nSta,:)=rho_jab_phasic(nS,restrictLagsvec);
    dec_bySub_Lag_mean(nS,nSta,:)=100*(rho_real_phasic(nS,restrictLagsvec)>rho_jab_phasic(nS,restrictLagsvec));
    %     for nSeg=1:length(begS)
    %         rho_att_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_att_Seg{nS}(trialphasic{nS},:,1),2)));
    %         rho_ign_bySub_Seg_mean(nS,nSta,nSeg)=mean(squeeze(mean(tem_ign_Seg{nS}(trialphasic{nS},:,2),2)));
    %         dec_bySub_Seg_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_Seg{nS}(trialtonic{nS},:,1),2)>mean(tem_ign_Seg{nS}(trialtonic{nS},:,1),2));
    %     end
end
