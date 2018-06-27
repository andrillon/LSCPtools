function peak_locs=detect_heartbeat(D,method,plotFlag)

if nargin<3
    plotFlag=0;
end
if nargin<2
    method=1;
end

%% pre-processing ECG
tm=D.time;
ECG_idx=match_str(D.chanlabels,'ECG');
ECG=squeeze(D(ECG_idx,:,1));
ECG_filt = ft_preproc_bandstopfilter(ECG,D.fsample,[45 55],5,'but','twopass');
ECG_filt=bandpass(ECG_filt,D.fsample,0.05,40,4);

if method==1
    ECG_norm=(ECG_filt-nanmedian(ECG_filt))./(nanstd(ECG_filt));
%% detect peak in ECG
    [pks,peak_locs] = findpeaks(ECG_norm);
    % clean small peaks
    % peak_locs(ECG_norm(peak_locs)<3)=[];
    % clean peaks too close
    limitHB=150;
    peak_locs(find(diff(peak_locs/D.fsample)<(1/(limitHB/60)))+1)=[];
    
    %% summing-up detection
    min_interval=prctile(diff(peak_locs)/D.fsample,1);
    max_interval=prctile(diff(peak_locs)/D.fsample,99);
    mean_interval=mean(diff(peak_locs)/D.fsample);
    fprintf('... ... %g heartbeats detected with an average HB rate of %2.1f hbm\n',length(peak_locs),1/mean_interval*60)
    fprintf('... ... intervals 1st percentile: %g / 99th percentile:  %g\n',min_interval,max_interval)
    
elseif method==2
    wt = modwt(ECG_filt,5);
    wtrec = zeros(size(wt));
    wtrec(4:5,:) = wt(4:5,:);
    y = imodwt(wtrec,'sym4');
    
    [pks,peak_locs] = findpeaks(y,'MinPeakHeight',5,...
        'MinPeakDistance',0.150);
    
    min_interval=prctile(diff(peak_locs)/D.fsample,1);
    max_interval=prctile(diff(peak_locs)/D.fsample,99);
    mean_interval=mean(diff(peak_locs)/D.fsample);
    fprintf('... ... %g heartbeats detected with an average HB rate of %2.1f hbm\n',length(peak_locs),1/mean_interval*60)
    fprintf('... ... intervals 1st percentile: %g / 99th percentile:  %g\n',min_interval,max_interval)
    
end
%%
if plotFlag
    figure;
    subplot(2,2,1:2); hold on; format_fig
    plot(ECG_norm)
    scatter(peak_locs,ECG_norm(peak_locs),'or')
    title('Heartbeat detection')
    
    % ERP on ECG
    subplot(2,2,3); hold on; format_fig
    ecg_peaks=nan(length(peak_locs),length((-0.1*D.fsample:0.4*D.fsample)));
    for k=1:length(peak_locs)
        ecg_peaks(k,:)=ECG_norm((-0.1*D.fsample:0.4*D.fsample)+peak_locs(k));
    end
    plot(-0.1:1/D.fsample:0.4,nanmean(ecg_peaks))
    title('ERP on ECG')
    xlim([-0.1 0.4])
    
    % ERP on Cz
    Cz_ecg_peaks=nan(length(peak_locs),length((-0.1*D.fsample:0.4*D.fsample)));
    fprintf('%3.0f%%\n',0)
    for k=1:length(peak_locs)
        fprintf('\b\b\b\b\b%3.0f%%\n',round(k/length(peak_locs)*100))
        Cz_ecg_peaks(k,:)=D(match_str(D.chanlabels,{'Cz'}),(-0.1*D.fsample:0.4*D.fsample)+peak_locs(k),1);
    end
    mytimes=-0.1:1/D.fsample:0.4;
    Cz_ecg_peaks=Cz_ecg_peaks-repmat(nanmean(Cz_ecg_peaks(:,mytimes<0),2),1,size(Cz_ecg_peaks,2));
    subplot(2,2,4); hold on; format_fig
    plot(mytimes,nanmean(Cz_ecg_peaks))
    title('ERP on Cz')
    xlim([-0.1 0.4])
end
