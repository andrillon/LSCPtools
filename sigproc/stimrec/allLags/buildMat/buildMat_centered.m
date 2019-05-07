rho_att_bySub_mean=[];rho_ign_bySub_mean=[];dec_bySub_mean=[];
rho_att_bySub_Lag=[];rho_ign_bySub_Lag=[];dec_bySub_Lag_mean=[];
rho_att_bySub_Seg_mean=[];rho_ign_bySub_Seg_mean=[];dec_bySub_Seg_mean=[];
goodtrial={};goodsubj={};
typeREMs={'burst','isolated','all'};
retrictLags=[0:25];restrictLagsvec=find(ismember(lagvec,retrictLags));
newlags=lagvec(restrictLagsvec);
lenS=length(begS);

for tpREM=1:length(typeREMs)
 count=0;
   typeREM=typeREMs{tpREM};
    for nS=find(cellfun(@(x) ~isempty(x),REM_compute))
        if strcmp(typeREM,'burst')
            restrictindex=REM_compute{nS}(:,5)==1;
        elseif strcmp(typeREM,'isolated')
            restrictindex=REM_compute{nS}(:,5)==0;
        elseif strcmp(typeREM,'all')
            restrictindex=ones(length(REM_compute{nS}),1);
        elseif strcmp(typeREM,'onsetBurst')
            restrictindex=REM_compute{nS}(:,5)==1 & REM_compute{nS}(:,6)==1;
        end
        if sum(restrictindex)>4
            count=count+1;
            goodtrial{tpREM}(count)=sum(restrictindex);
            goodsubj{tpREM}(count)=nS;
        else
            continue
        end
        rho_att_bySub{nS}=squeeze(mean(tem_att{nS}(restrictindex,restrictLagsvec,:),2));
        rho_ign_bySub{nS}=squeeze(mean(tem_ign{nS}(restrictindex,restrictLagsvec,:),2));
        dec_bySub{nS}=100*mean(rho_att_bySub{nS}>rho_ign_bySub{nS},1);
        rho_att_Lag{nS}=squeeze(mean(tem_att{nS}(restrictindex,restrictLagsvec,:),3));
        rho_ign_Lag{nS}=squeeze(mean(tem_ign{nS}(restrictindex,restrictLagsvec,:),3));
        dec_Lag{nS}=100*mean(rho_att_Lag{nS}>rho_ign_Lag{nS},1);

    end
    rho_att_bySub_Seg_mean(goodsubj{tpREM},tpREM,:)=reshape(cell2mat(cellfun(@(x) mean(x,1),rho_att_bySub(goodsubj{tpREM}),'uni',0)),length(goodsubj{tpREM}),lenS);
    rho_ign_bySub_Seg_mean(goodsubj{tpREM},tpREM,:)=reshape(cell2mat(cellfun(@(x) mean(x,1),rho_ign_bySub(goodsubj{tpREM}),'uni',0)),length(goodsubj{tpREM}),lenS);
    dec_bySub_Seg_mean(goodsubj{tpREM},tpREM,:)=reshape(cell2mat(cellfun(@(x) mean(x,1),dec_bySub(goodsubj{tpREM}),'uni',0)),length(goodsubj{tpREM}),lenS);
    
    rho_att_bySub_mean(goodsubj{tpREM},tpREM)=mean(rho_att_bySub_Seg_mean(goodsubj{tpREM},tpREM),2);
    rho_ign_bySub_mean(goodsubj{tpREM},tpREM)=mean(rho_ign_bySub_Seg_mean(goodsubj{tpREM},tpREM),2);
    dec_bySub_mean(goodsubj{tpREM},tpREM)=mean(dec_bySub_Seg_mean(goodsubj{tpREM},tpREM),2);
    
    rho_att_bySub_Lag(goodsubj{tpREM},tpREM,1:length(restrictLagsvec))=reshape(cell2mat(cellfun(@(x) mean(x,1),rho_att_Lag(goodsubj{tpREM}),'uni',0)),length(goodsubj{tpREM}),length(restrictLagsvec));
    rho_ign_bySub_Lag(goodsubj{tpREM},tpREM,1:length(restrictLagsvec))=reshape(cell2mat(cellfun(@(x) mean(x,1),rho_ign_Lag(goodsubj{tpREM}),'uni',0)),length(goodsubj{tpREM}),length(restrictLagsvec));
    dec_bySub_Lag_mean(goodsubj{tpREM},tpREM,1:length(restrictLagsvec))=reshape(cell2mat(cellfun(@(x) mean(x,1),dec_Lag(goodsubj{tpREM}),'uni',0)),length(goodsubj{tpREM}),length(restrictLagsvec));
    
end
NameStages=typeREMs;
