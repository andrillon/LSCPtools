goodsubj=cell(1,length(NameStages));
Rmatrix_diff_att_ign=[];Rmatrix_diff_att_ign_Seg=[];
baselineERP=[-0.2,0];
window_ERP_plot=[-0.5,1];
baseline_idx=find(ismember(window_ERP(1)*fs:window_ERP(2)*fs,baselineERP(1)*fs:baselineERP(2)*fs));
plot_idx=find(ismember(window_ERP(1)*fs:window_ERP(2)*fs,window_ERP_plot(1)*fs:window_ERP_plot(2)*fs));
att_meanERP_bySub=[];ign_meanERP_bySub=[];
att_meanERP_bySub_half=[];ign_meanERP_bySub_half=[];
ttestplot=1;

for nSta=1:length(NameStages)
    for nS=1:length(subject_id)
        if count(nS,nSta)>=4
            att_baselinecorr=att_meanERP{nS,nSta}(:,plot_idx)-repmat(mean(att_meanERP{nS,nSta}(:,baseline_idx),2),1,length(plot_idx));
            ign_baselinecorr=ign_meanERP{nS,nSta}(:,plot_idx)-repmat(mean(ign_meanERP{nS,nSta}(:,baseline_idx),2),1,length(plot_idx));
            att_meanERP_bySub(nS,nSta,:)=nanmean(att_baselinecorr,1);
            ign_meanERP_bySub(nS,nSta,:)=nanmean(ign_baselinecorr,1);
                        
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
