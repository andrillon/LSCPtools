NameStages={'tREM','pREM'};
clear goodsubj
goodsubj{1}=find(lenphases(:,1)>5*60);%minimal 5 minutes
goodsubj{2}=find(lenphases(:,2)>5*60);%minimal 5 minutes

rho_att_bySub_mean=zeros(length(subject_id),2);rho_ign_bySub_mean=zeros(length(subject_id),2);dec_bySub_mean=zeros(length(subject_id),2);
rho_att_bySub_Lag=zeros(length(subject_id),2,length(lagvec));rho_ign_bySub_Lag=zeros(length(subject_id),2,length(lagvec));dec_bySub_Lag_mean=zeros(length(subject_id),2,length(lagvec));
trialphasic=cell(1,length(subject_id));trialtonic=cell(1,length(subject_id));
begS=[];begS_half=[];

for nS=1:length(subject_id)
    
    
    %% tonic
    nSta=1;
    rho_att_bySub_mean(nS,nSta)=rho_real_tonic(nS);
    rho_ign_bySub_mean(nS,nSta)=rho_jab_tonic(nS);
    dec_bySub_mean(nS,nSta)=100*(rho_att_bySub_mean(nS,nSta)>rho_ign_bySub_mean(nS,nSta));
    
    rho_att_bySub_Lag(nS,nSta,:)=rho_real_tonic(nS);
    rho_ign_bySub_Lag(nS,nSta,:)=rho_jab_tonic(nS);
    dec_bySub_Lag_mean(nS,nSta,:)=100*(rho_real_tonic(nS)>rho_jab_tonic(nS));
    
    %% phasic
    nSta=2;
    rho_att_bySub_mean(nS,nSta)=mean(rho_real_phasic(nS,restrictLagsvec));
    rho_ign_bySub_mean(nS,nSta)=mean(rho_jab_phasic(nS,restrictLagsvec));
    dec_bySub_mean(nS,nSta)=100*(rho_att_bySub_mean(nS,nSta)>rho_ign_bySub_mean(nS,nSta));
    
    rho_att_bySub_Lag(nS,nSta,:)=rho_real_phasic(nS);
    rho_ign_bySub_Lag(nS,nSta,:)=rho_jab_phasic(nS);
    dec_bySub_Lag_mean(nS,nSta,:)=100*(rho_real_phasic(nS)>rho_jab_phasic(nS));
end
