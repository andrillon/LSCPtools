goodsubj=cell(1,length(NameStages));
Rmatrix_diff_att_ign=[];Rmatrix_diff_dec=[];Rmatrix_diff_att_ign_Seg=[];Rmatrix_diff_dec_Seg=[];
rho_att_bySub_mean=[];rho_ign_bySub_mean=[];dec_bySub_mean=[];
rho_att_bySub_Lag=[];rho_ign_bySub_Lag=[];dec_bySub_Lag_mean=[];
rho_att_bySub_half_mean=[];rho_ign_bySub_half_mean=[];dec_bySub_half_mean=[];
rho_att_bySub_Seg_mean=[];rho_ign_bySub_Seg_mean=[];dec_bySub_Seg_mean=[];
ttestplot=1;

for nSta=1:length(NameStages)
    for nS=1:length(subject_id)
        if count(nS,nSta)>=4
            rho_att_bySub_mean(nS,nSta)=mean(mean(tem_att{nS,nSta},2));
            rho_ign_bySub_mean(nS,nSta)=mean(mean(tem_ign{nS,nSta},2));
            dec_bySub_mean(nS,nSta,:)=100*mean(mean(tem_att{nS,nSta},2)>mean(tem_ign{nS,nSta},2));
            
            rho_att_bySub_Lag(nS,nSta,:)=mean(tem_att{nS,nSta},1);
            rho_ign_bySub_Lag(nS,nSta,:)=mean(tem_ign{nS,nSta},1);
            dec_bySub_Lag_mean(nS,nSta,:)=100*mean(tem_att{nS,nSta}>tem_ign{nS,nSta},1);
            
            %             vecone=ones(count(nS,nSta),1);
            %             attrec=[mean(tem_att{nS,nSta},2),nSta*vecone,vecone,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
            %             ignrec=[mean(tem_att{nS,nSta},2),nSta*vecone,vecone*0,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
            %             decvec=[100*(mean(tem_att{nS,nSta},2)>mean(tem_ign{nS,nSta},2)),nSta*vecone,vecone,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
            %
            %             Rmatrix_diff_att_ign=[Rmatrix_diff_att_ign;attrec;ignrec];
            %             Rmatrix_diff_dec=[Rmatrix_diff_dec;decvec];
            if ~isempty(begS_half)
                for nSeg=1:length(begS_half)
                    rho_att_bySub_half_mean(nS,nSta,nSeg)=mean(mean(tem_att_half{nS,nSta}(:,:,nSeg),2));
                    rho_ign_bySub_half_mean(nS,nSta,nSeg)=mean(mean(tem_ign_half{nS,nSta}(:,:,nSeg),2));
                    dec_bySub_half_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_half{nS,nSta}(:,:,nSeg),2)>mean(tem_ign_half{nS,nSta}(:,:,nSeg),2));
                    
                    %                     attrecseg=[mean(tem_att_half{nS,nSta}(:,:,nSeg),2),nSta*vecone,vecone,vecone*nSeg,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
                    %                     ignrecseg=[mean(tem_ign_half{nS,nSta}(:,:,nSeg),2),nSta*vecone,vecone*0,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
                    %                     decvec_half=[100*(mean(tem_att_half{nS,nSta}(:,:,nSeg),2)>mean(tem_ign_half{nS,nSta}(:,:,nSeg),2)),nSta*vecone,vecone,vecone*nSeg,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
                    
                    %                     Rmatrix_diff_att_ign_half=[Rmatrix_diff_att_ign_half;attrecseg;ignrecseg];
                    %                     Rmatrix_diff_dec_half=[Rmatrix_diff_dec_half;decvec_half];
                end
            end
            
            if ~isempty(begS)
                for nSeg=1:length(begS)
                    rho_att_bySub_Seg_mean(nS,nSta,nSeg)=mean(mean(tem_att_Seg{nS,nSta}(:,:,nSeg),2));
                    rho_ign_bySub_Seg_mean(nS,nSta,nSeg)=mean(mean(tem_ign_Seg{nS,nSta}(:,:,nSeg),2));
                    dec_bySub_Seg_mean(nS,nSta,nSeg)=100*mean(mean(tem_att_Seg{nS,nSta}(:,:,nSeg),2)>mean(tem_ign_Seg{nS,nSta}(:,:,nSeg),2));
                    
                    %                     attrecseg=[mean(tem_att_Seg{nS,nSta}(:,:,nSeg),2),nSta*vecone,vecone,vecone*nSeg,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
                    %                     ignrecseg=[mean(tem_ign_Seg{nS,nSta}(:,:,nSeg),2),nSta*vecone,vecone*0,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
                    %                     decvec_Seg=[100*(mean(tem_att_Seg{nS,nSta}(:,:,nSeg),2)>mean(tem_ign_Seg{nS,nSta}(:,:,nSeg),2)),nSta*vecone,vecone,vecone*nSeg,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
                    
                    %                     Rmatrix_diff_att_ign_Seg=[Rmatrix_diff_att_ign_Seg;attrecseg;ignrecseg];
                    %                     Rmatrix_diff_dec_Seg=[Rmatrix_diff_dec_Seg;decvec_Seg];
                end
            end
        end
    end
    goodsubj{nSta}=find(count(:,nSta)>4);
end
