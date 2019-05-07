%% plot Results
begS=[];contmax=60;fltrate=0.2;
baselineTime=[-3,0];%baseline in seconds
NameStages={'RealTr','JabbTr','Wake','N2','REM'};
goodsubj=cell(1,length(NameStages));
Rmat_Zv_trial=[];Rmat_Zv_baseline=[];Rmatrix_diff_dec=[];
Zv_mean=zeros(length(subject_id),length(NameStages));
Zv_baseline_mean=zeros(length(subject_id),length(NameStages));
Zv_Seg_mean=zeros(length(subject_id),length(NameStages),length(begS_half));
Zv_mean_cont=zeros(length(subject_id),length(NameStages),contmax/fltrate);
for nSta=1:length(NameStages)
    for nS=1:length(subject_id)
        %             mean across trials
        if count(nS,nSta)>0
            Zv_mean(nS,nSta)=mean(Lempel_Ziv_trial_Zv{nS,nSta}(:,1)./Lempel_Ziv_trial_perm{nS,nSta}(:,1));
            Zv_baseline_mean(nS,nSta)=mean(Lempel_Ziv_baseline_Zv{nS,nSta}(:,1)./Lempel_Ziv_baseline_perm{nS,nSta}(:,1));
            Zv_mean_cont(nS,nSta,:)=squeeze(mean(Lempel_Ziv_trial_norm_continuous{nS,nSta}(:,1,1:contmax/fltrate),1));
%             vecone=ones(count(nS,nSta),1);
% %             Zv_trial=[Lempel_Ziv_trial_norm{nS,nSta}(:,1),nSta*vecone,vecone,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
% %             Zv_baseline=[Lempel_Ziv_baseline_norm{nS,nSta}(:,1),nSta*vecone,vecone*0,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
% %             Rmat_Zv_trial=[Rmat_Zv_trial;Zv_trial];Rmat_Zv_baseline=[Rmat_Zv_baseline;Zv_baseline];
%             
            if ~isempty(begS_half)
                for nSeg=1:length(begS_half)
                    Zv_Seg_mean(nS,nSta,nSeg)=mean(Lempel_Ziv_trial_norm_nSeg{nS,nSta}(:,1,nSeg));
%                    Zv_trial_Seg=[Lempel_Ziv_trial_norm_nSeg{nS,nSta}(:,nSeg),nSeg*vecone,nSta*vecone,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
                end
            end
        end

        goodsubj{nSta}=find(count(:,nSta)>=4);
    end
end