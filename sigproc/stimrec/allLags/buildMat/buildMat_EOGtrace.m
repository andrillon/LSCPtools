lEOG_bySub_mean={};rEOG_bySub_mean=[];dec_bySub_mean=[];
lEOG_bySub_Lag={};rEOG_bySub_Lag=[];dec_bySub_Lag_mean=[];
lEOG_bySub_Seg_mean=[];rEOG_bySub_Seg_mean=[];dec_bySub_Seg_mean=[];
goodtrial={};goodsubj={};
typeREMs={'onsetBurst','burst','isolated'};

lenS=length(begS);

for tpREM=1:length(typeREMs)
 count=0;
   typeREM=typeREMs{tpREM};
    for nS=1:length(subject_id)
        try
        if strcmp(typeREM,'burst')
            restrictindex=REM_compute{nS}(restrict_trials_EOG{nS},5)==1
        elseif strcmp(typeREM,'isolated')
            restrictindex=REM_compute{nS}(restrict_trials_EOG{nS},5)==0;
        elseif strcmp(typeREM,'all')
            restrictindex=ones(length(REM_compute{nS}),1);
        elseif strcmp(typeREM,'onsetBurst')
            restrictindex=REM_compute{nS}(restrict_trials_EOG{nS},5)==1 & REM_compute{nS}(restrict_trials_EOG{nS},6)==1;
        end
        catch
            restrictindex=[];
        end
        if sum(restrictindex)>4
            count=count+1;
            goodtrial{tpREM}(count)=sum(restrictindex);
            goodsubj{tpREM}(count)=nS;
        else
            continue
        end
        lEOGfilt=squeeze(lEOG_trace{nS}(restrictindex,1,:));
        rEOGfilt=squeeze(rEOG_trace{nS}(restrictindex,1,:));
        
%         filt1=ft_preproc_highpassfilter(squeeze(lEOG_trace{nS}(restrictindex,1,:)),500,filterEOG(1),5,'but');%,'onepass-zerophase');
%         lEOGfilt=ft_preproc_lowpassfilter(filt1,500,filterEOG(2),5,'but');%,'onepass-zerophase');
%         filt1=ft_preproc_highpassfilter(squeeze(rEOG_trace{nS}(restrictindex,1,:)),500,filterEOG(1),5,'but');%,'firws','onepass-zerophase');
%         rEOGfilt=ft_preproc_lowpassfilter(filt1,500,filterEOG(2),5,'but');%,'onepass-zerophase');
% %   
        lEOG_bySub{nS}=squeeze(mean(lEOGfilt,1));
        rEOG_bySub{nS}=squeeze(mean(rEOGfilt,1));

    end    
    lEOG_bySub_mean{tpREM}=cell2mat(cellfun(@(x) mean(x,1),lEOG_bySub(goodsubj{tpREM}),'uni',0)');
    rEOG_bySub_mean{tpREM}=cell2mat(cellfun(@(x)  mean(x,1),rEOG_bySub(goodsubj{tpREM}),'uni',0)');
    
end
NameStages=typeREMs;
