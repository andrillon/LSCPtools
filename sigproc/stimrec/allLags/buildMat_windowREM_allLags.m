NameStages={'tREM','pREM'};
rho_att_bySub_mean=zeros(length(subject_id),1);rho_ign_bySub_mean=zeros(length(subject_id),1);dec_bySub_mean=zeros(length(subject_id),1);
trialphasic=cell(1,length(subject_id));trialtonic=cell(1,length(subject_id));
rho_att_bySub_Lag=zeros(length(subject_id),length(lagvec));rho_ign_bySub_Lag=zeros(length(subject_id),length(lagvec));dec_bySub_Lag_mean=zeros(length(subject_id),length(lagvec));
goodsubj=cell(1,2);
for nS=1:length(subject_id)
    trialphasic{nS}=find(phaseREM{nS}==6);
    trialtonic{nS}=find(phaseREM{nS}==5);
    if length(trialtonic{nS})>4
        %% tonic
        nSta=1;
        rho_att_bySub_mean(nS,nSta)=mean(tem_att{nS}(trialtonic{nS}));
        rho_ign_bySub_mean(nS,nSta)=mean(tem_ign{nS}(trialtonic{nS}));
        dec_bySub_mean(nS,nSta)=100*mean(tem_att{nS}(trialtonic{nS})>tem_ign{nS}(trialtonic{nS}));
    end
    %% phasic
    if length(trialphasic{nS})>4
        nSta=2;
        rho_att_bySub_mean(nS,nSta)=mean(tem_att{nS}(trialphasic{nS}));
        rho_ign_bySub_mean(nS,nSta)=mean(tem_ign{nS}(trialphasic{nS}));
        dec_bySub_mean(nS,nSta)=100*mean(tem_att{nS}(trialphasic{nS})>mean(tem_ign{nS}(trialphasic{nS})));
    end
end
goodsubj{1}=find(cellfun(@length,trialtonic)>=4);
goodsubj{2}=find(cellfun(@length,trialphasic)>=4);
