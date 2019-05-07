%% plot Results
load([savePath 'trialREM',num2str(filterband(1)),'to',num2str(filterband(2)),'.mat'],'nE_id','nrREMtrial','subject_id')
nE_id_REM=nE_id;subject_id_REM=subject_id;
load([savePath 'LempelZiv_result.mat'],'nE_id','subject_id')
nE_id_Lpz=nE_id;
NameStages={'tREM','pREM'};
goodsubj=cell(1,length(NameStages));
Rmat_Zv_trial=[];Rmat_Zv_baseline=[];Rmatrix_diff_dec=[];
Zv_mean=zeros(length(subject_id_REM),length(NameStages));
Zv_baseline_mean=zeros(length(subject_id_REM),length(NameStages));
Zv_Seg_mean=zeros(length(subject_id_REM),length(NameStages),length(begS_half));
Zv_mean_cont=zeros(length(subject_id_REM),length(NameStages),contmax/fltrate);

for nS=1:length(subject_id_REM)
    nS_Lpz=find(cellfun(@(x) strcmp(x,subject_id_REM{nS}),subject_id));
    if ~isempty(nS_Lpz)
        [idfind,idREM]=ismember(nE_id_Lpz{nS_Lpz,end},nE_id_REM{nS});
        REMphasictrial{nS}=find(nrREMtrial{nS}(idREM)>0);
        REMtonictrial{nS}=find(nrREMtrial{nS}(idREM)==0);
        
        if ~isempty(REMtonictrial{nS})
            Zv_mean(nS,1)=mean(Lempel_Ziv_trial_Zv{nS_Lpz,end}(REMtonictrial{nS},1)./Lempel_Ziv_trial_perm{nS_Lpz,end}(REMtonictrial{nS},1));
            Zv_baseline_mean(nS,1)=mean(Lempel_Ziv_baseline_Zv{nS_Lpz,end}(REMtonictrial{nS},1)./Lempel_Ziv_baseline_perm{nS_Lpz,end}(REMtonictrial{nS},1));
            Zv_mean_cont(nS,1,:)=squeeze(mean(Lempel_Ziv_trial_norm_continuous{nS_Lpz,end}(REMtonictrial{nS},1,1:contmax/fltrate),1));
        end
        
        if ~isempty(REMphasictrial{nS})
            Zv_mean(nS,2)=mean(Lempel_Ziv_trial_Zv{nS_Lpz,end}(REMphasictrial{nS},1)./Lempel_Ziv_trial_perm{nS_Lpz,end}(REMphasictrial{nS},1));
            Zv_baseline_mean(nS,2)=mean(Lempel_Ziv_baseline_Zv{nS_Lpz,end}(REMphasictrial{nS},1)./Lempel_Ziv_baseline_perm{nS_Lpz,end}(REMphasictrial{nS},1));
            Zv_mean_cont(nS,2,:)=squeeze(mean(Lempel_Ziv_trial_norm_continuous{nS_Lpz,end}(REMphasictrial{nS},1,1:contmax/fltrate),1));
        end
        %             vecone=ones(count(nS,nSta),1);
        % %             Zv_trial=[Lempel_Ziv_trial_norm{nS,nSta}(:,1),nSta*vecone,vecone,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
        % %             Zv_baseline=[Lempel_Ziv_baseline_norm{nS,nSta}(:,1),nSta*vecone,vecone*0,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
        % %             Rmat_Zv_trial=[Rmat_Zv_trial;Zv_trial];Rmat_Zv_baseline=[Rmat_Zv_baseline;Zv_baseline];
        %
        if ~isempty(begS_half)
            for nSeg=1:length(begS_half)
                if ~isempty(REMtonictrial{nS})
                    Zv_Seg_mean(nS,1,nSeg)=mean(Lempel_Ziv_trial_norm_nSeg{nS_Lpz,end}(REMtonictrial{nS},1,nSeg));
                end
                if ~isempty(REMphasictrial{nS})
                    Zv_Seg_mean(nS,2,nSeg)=mean(Lempel_Ziv_trial_norm_nSeg{nS_Lpz,end}(REMphasictrial{nS},2,nSeg));
                end
                %                    Zv_trial_Seg=[Lempel_Ziv_trial_norm_nSeg{nS,nSta}(:,nSeg),nSeg*vecone,nSta*vecone,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
            end
        end
    end
    
end
% goodsubj={find(cellfun(@isempty,REMtonictrial)==0),find(cellfun(@isempty,REMphasictrial)==0)};
goodsubj={find(cellfun(@length,REMtonictrial)>=4),find(cellfun(@length,REMphasictrial)>=4)};
