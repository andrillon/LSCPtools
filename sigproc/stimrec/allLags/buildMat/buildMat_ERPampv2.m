goodsubj=cell(1,length(NameStages));
Rmatrix_diff_att_ign=[];Rmatrix_diff_att_ign_Seg=[];

rho_att_bySub_mean=[];rho_ign_bySub_mean=[];
rho_att_bySub_Lag=[];rho_ign_bySub_Lag=[];dec_bySub_Lag_mean=[];
rho_att_bySub_half_mean=[];rho_ign_bySub_half_mean=[];dec_bySub_half_mean=[];
rho_att_bySub_Seg_mean=[];rho_ign_bySub_Seg_mean=[];dec_bySub_Seg_mean=[];
ttestplot=1;
baselinewindow=[-0.5,0];
ERPwindow=[-0.5,1];
begS=[0,30];endS=[30,60];
for nSta=1:length(NameStages)
    for nS=1:length(subject_id)
        if count(nS,nSta)>=4
            
            %filtering
            filtering_att_low=ft_preproc_highpassfilter(att_meanERP{nS,nSta},500,filterERP(1),[],'firws','onepass-zerophase');
            filtering_att_filt=ft_preproc_lowpassfilter(filtering_att_low,500,filterERP(2),[],'firws','onepass-zerophase');
            filtering_ign_low=ft_preproc_highpassfilter(att_meanERP{nS,nSta},500,filterERP(1),[],'firws','onepass-zerophase');
            filtering_ign_filt=ft_preproc_lowpassfilter(filtering_ign_low,500,filterERP(2),[],'firws','onepass-zerophase');
            
            %ERP and baseline correction
            baselineidx=(baselinewindow(1)-window_lowpass(1))*fs:(baselinewindow(2)-window_lowpass(1))*fs;
            ERPidx=(ERPwindow(1)-window_lowpass(1))*fs:(ERPwindow(2)-window_lowpass(1))*fs;
            att_ERP_bySub=filtering_att_filt(ERPidx,:)-repmat(mean(filtering_att_high(baselineidx,:),2),1,length(ERPidx));
            att_meanERP_bySub(nS,nSta,:)=nanmean(att_ERP_bySub,1);
            ign_ERP_bySub=filtering_ign_filt(ERPidx,:)-repmat(mean(filtering_ign_filt(baselineidx,:),2),1,length(ERPidx));
            ign_meanERP_bySub(nS,nSta,:)=nanmean(ign_ERP_bySub,1);
                        
            if ~isempty(begS_half)
                for nSeg=1:length(begS_half)
                    Segix=[AvSilenceNumber{nS,nSta,1}]>begS(nSeg) & [AvSilenceNumber{nS,nSta,1}]<endS(nSeg);
                    att_meanERP_bySub_half(nS,nSta,nSeg,:)=nanmean(att_ERP_bySub(Segix,:),1);
                    ign_meanERP_bySub_half(nS,nSta,nSeg,:)=nanmean(ign_ERP_bySub(Segix,:),1);
                end
            end
            
        end
    end
    goodsubj{nSta}=find(count(:,nSta)>4);
end
