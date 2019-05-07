%% parameters for reconstruction
Fs=500;%frequency sampling of EEG to lock the trial 
fltrate=1/Fs;
begS=[1,30./fltrate];endS=[30/fltrate+1,60/fltrate];%time windows for 0-30s,30s-60s
cont_max=60;%maximum time in seconds for continuous 
baseline=[-3,0];%baseline in seconds
NameStages={'RealTr','JabbTr','W','N2','REM'};
maxdiff=[];
%time windows for 0-30s,30s-60s
begS_half=[1,30./fltrate];endS_half=[30/fltrate+1,60/fltrate];

%initialisation
count=zeros(length(subject_id),length(NameStages));
Lempel_Ziv_trial_Zv=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_perm=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_norm=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_HL=cell(length(subject_id),length(NameStages));
Lempel_Ziv_baseline_Zv=cell(length(subject_id),length(NameStages));
Lempel_Ziv_baseline_perm=cell(length(subject_id),length(NameStages));
Lempel_Ziv_baseline_norm=cell(length(subject_id),length(NameStages));
Lempel_Ziv_baseline_HL=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_Zv_nSeg=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_perm_nSeg=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_norm_nSeg=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_HL_nSeg=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_Zv_continuous=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_perm_continuous=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_norm_continuous=cell(length(subject_id),length(NameStages));
Lempel_Ziv_trial_HL_continuous=cell(length(subject_id),length(NameStages));
baselineTime=[-4,0];
problem_nE=[];

%load Lempel Ziv data
load([ExtractedDataPath,'LempelZiv2.mat']);

windowsize=LempelZiv(1).windowsize;
%% 
for nS=17:length(subject_id)
    
    %% subject def
    subName=[Expe subject_id{nS}];
    SubID=subject_id{nS};
    fprintf('... Subject: %s\n',SubID);
    if ~ismember(SubID,{LempelZiv.subjID})
        continue
    end
    %% Retrive behavioral data
    SubjectBehavData=load([BehavPath filesep subName filesep 'Result_' subName '.mat']);
    fprintf('... ... behavioral data retrived\n')
    
    %% find lempel ziv data of the subject
    id_nS=find(cellfun(@(x) strcmp(x,SubID),{LempelZiv.subjID}));
    begS_time=LempelZiv(nS).timewindow(1,:);
    endS_time=LempelZiv(nS).timewindow(2,:);
    %% Retrieve events EEG data
    try
        load([ExtractedDataPath,'data_extracted' SubID '.mat'],'myevents')
    catch
        myevents=ft_read_event([rawEEGPath,subName,'.raw']);
    end
    
    %% retrieve structure of data : ordered by type of events, with scoring in tagSleep
    
    CPR_trialevents
    CPR_trialstructure
    trial_sleep_scoring

    %% 
    for st=1:4
        for nT=1:length([trials_st(st).id])
            if st==4
                if trials_st(st).tagSleepScoring(nT)==2;
                    nSta=4;
                elseif trials_st(st).tagSleepScoring(nT)>=5
                    nSta=5;
                else
                    continue
                end
            else
                nSta=st;
            end
            
            nE=trials_st(st).id(nT);
            fprintf('.. %g/40 ..\n',nE)
            count(nS,nSta)=count(nS,nSta)+1;
            nE_id{nS,nSta}(count(nS,nSta))=nE;
            
            window_time=[trials_st(st).trials(nT,1),trials_st(st).trials(nT,2)];
            timewindow_trial=find(begS_time>window_time(1) & endS_time<window_time(2));
            timewindow_baseline=find(begS_time>window_time(1)+baselineTime(1)*500 & endS_time<window_time(1)+baselineTime(2)*500);
            
            %trial
            Lempel_Ziv_trial_Zv{nS,nSta}(count(nS,nSta),1:3)=mean(LempelZiv(id_nS).Zv(timewindow_trial,:),1);
            Lempel_Ziv_trial_perm{nS,nSta}(count(nS,nSta),1:3)=mean(LempelZiv(id_nS).perm(timewindow_trial,:),1);
            Lempel_Ziv_trial_norm{nS,nSta}(count(nS,nSta),1:3)=mean(LempelZiv(id_nS).Zv(timewindow_trial,1)./LempelZiv(id_nS).perm(timewindow_trial,1),1);
            Lempel_Ziv_trial_HL{nS,nSta}(count(nS,nSta))=mean(LempelZiv(id_nS).HL(timewindow_trial));
            
            %baseline
            Lempel_Ziv_baseline_Zv{nS,nSta}(count(nS,nSta),:)=mean(LempelZiv(id_nS).Zv(timewindow_baseline,:),1);
            Lempel_Ziv_baseline_perm{nS,nSta}(count(nS,nSta),:)=mean(LempelZiv(id_nS).perm(timewindow_baseline,:),1);
            Lempel_Ziv_baseline_norm{nS,nSta}(count(nS,nSta),:)=mean(LempelZiv(id_nS).Zv(timewindow_baseline,1)./LempelZiv(id_nS).perm(timewindow_baseline,1),1);
            Lempel_Ziv_baseline_HL{nS,nSta}(count(nS,nSta))=mean(LempelZiv(id_nS).HL(timewindow_baseline));
            
            %30s-60s
            if ~isempty(begS_half)
                for nSeg=1:length(begS_half)
                    timewindow_trial=find(begS_time>window_time(1)+begS(nSeg) & endS_time<window_time(1)+endS(nSeg));
                    Lempel_Ziv_trial_Zv_nSeg{nS,nSta}(count(nS,nSta),1:3,nSeg)=mean(LempelZiv(id_nS).Zv(timewindow_trial,:),1);
                    Lempel_Ziv_trial_perm_nSeg{nS,nSta}(count(nS,nSta),1:3,nSeg)=mean(LempelZiv(id_nS).perm(timewindow_trial,:),1);
                    Lempel_Ziv_trial_norm_nSeg{nS,nSta}(count(nS,nSta),1:3,nSeg)=mean(LempelZiv(id_nS).Zv(timewindow_trial,1)./LempelZiv(id_nS).perm(timewindow_trial,1),1);
                    Lempel_Ziv_trial_HL_nSeg{nS,nSta}(count(nS,nSta),nSeg)=mean(LempelZiv(id_nS).HL(timewindow_trial));
                end
            end
            
            %continuous
            timewindow_trial=find(begS_time>window_time(1)+baselineTime(1)*500 & endS_time<window_time(1)+cont_max*Fs);
            for nSeg=1:length(timewindow_trial)
                Lempel_Ziv_trial_Zv_continuous{nS,nSta}(count(nS,nSta),1:3,nSeg)=LempelZiv(id_nS).Zv(timewindow_trial(nSeg));
                Lempel_Ziv_trial_perm_continuous{nS,nSta}(count(nS,nSta),1:3,nSeg)=LempelZiv(id_nS).perm(timewindow_trial(nSeg));
                Lempel_Ziv_trial_norm_continuous{nS,nSta}(count(nS,nSta),1:3,nSeg)=LempelZiv(id_nS).Zv(timewindow_trial(nSeg))./LempelZiv(id_nS).perm(timewindow_trial(nSeg));
                Lempel_Ziv_trial_HL_continuous{nS,nSta}(count(nS,nSta),nSeg)=LempelZiv(id_nS).HL(timewindow_trial(nSeg));
            end
        end
    end
end
