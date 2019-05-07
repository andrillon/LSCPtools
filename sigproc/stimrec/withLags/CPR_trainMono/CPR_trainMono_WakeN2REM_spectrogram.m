%% initialisation
subject_id = {'006','008','009','010','011','012','013','014',...
    '015','016','017','019','020','021','022','023','024','026',...
    '027','028','029','030','031','032','033','034','035','036',...
    '037','038','039'};%'011',;%,'050'};%,'090','100','130','150','160','170','180','200','210','220'};%,'190'};

Expe = 'CPR';
userID='Matthieu';

%%% path define
if strcmp(userID,'Celia')
    CPR_path_config_Celia
elseif strcmp(userID,'Matthieu')
    pathGit =genpath('/home/lscp00/Work/EEGgit/Sommeil/')
    addpath(pathGit)
    CPR_path_define_Matthieu
elseif strcmp(userID,'Thomas')
    CPR_path_define_Thomas
end


%%% parameters for reconstruction
fltrate=0.01;%sampled at 10Hz
begS=[1,30./fltrate];endS=[30/fltrate+1,60/fltrate];%time windows for 0-30s,30s-60s
Fs=500;

NameStages={'W'};%,'N2','REM'};
Numstages={[0,1]};%,2,[5,6]};
maxdiff=[];
count=zeros(length(subject_id),length(Numstages));
tem_att=cell(length(subject_id),length(Numstages));
tem_ign=cell(length(subject_id),length(Numstages));

chanSelec=[1:65];
bandfreq=[0.5,2;4,8;8,12;12,18];
side_Tale=cell(length(subject_id),length(Numstages));
nE_id=cell(length(subject_id),length(Numstages));
problem_nE=[];

%% Reconstruction
%10,18,24,29
%2-8:25,27,29

for nS=1:length(subject_id)
    
    %%% Import in SPM
    subName=[Expe subject_id{nS}];
    SubID=subject_id{nS};
    fprintf('... Subject: %s\n',SubID);
    
    %%% Retrive behavioral data
    SubjectBehavData=load([BehavPath filesep subName filesep 'Result_' subName '.mat']);
    fprintf('... ... behavioral data retrived\n')
    
    %%% Retrieve events EEG data
    try
        load([ExtractedDataPath,'data_extracted' SubID '.mat'],'myevents','data_full')
    catch
        myevents=ft_read_event([rawEEGPath,subName,'.raw']);
        data_full=ft_read_data([rawEEGPath,subName,'.raw']);
        save([ExtractedDataPath,'data_extracted' SubID '.mat'],'myevents','data_full')
    end
    
    %%% retrieve structure of data : ordered by type of events, with
    %%% scoring in tagSleep
    CPR_trialevents
    CPR_trialstructure
    trial_sleep_scoring
    
    clear g
    %%% training the model
    fprintf('... ... training model on di-otic\n')
    train_data=[];train_stim=[];
    
    fprintf('... ... test model on di-chotic\n')
    nSta=1;
    %%
    %Forced Wake Trials
    for nT=1:length([trials_st(3).id])
        try
            %       fprintf('.. %g/40 ..\n',nE)
            nE=trials_st(3).id(nT);
            count(nS,nSta)=count(nS,nSta)+1;
            window_time=trials_st(3).trials(nT,1):trials_st(3).trials(nT,2);
            
            for bandf=1:length(bandfreq)
                freqrange1=bandfreq(bandf,:);
                freq_stats{nS,nSta}(count(nS,nSta),bandf,chanSelec) = bandpower(data_full(chanSelec,window_time)',Fs,freqrange1);
            end
            
            nE_id{nS,nSta}(count(nS,nSta))=nE;
            if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                side_Tale{nS,nSta}(count(nS,nSta))=1;
            elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                side_Tale{nS,nSta}(count(nS,nSta))=2;
            end
            
            %30s-60s
            if ~isempty(begS)
                for nSeg=1:length(begS)
                    window_time=[trials_st(3).trials(nT,1)+(nSeg-1)*30*Fs]:min([trials_st(3).trials(nT,1)+(nSeg)*30*Fs,trials_st(3).trials(nT,2)]);
                    
                    for bandf=1:length(bandfreq)
                        freqrange1=bandfreq(bandf,:);
                        freq_stats_Seg{nS,nSta}(count(nS,nSta),nSeg,bandf,chanSelec) = bandpower(data_full(chanSelec,window_time)',Fs,bandfreq(bandf,:));
                    end
                end
            end
        catch
            problem_nE=[problem_nE;nS,nE]
        end
        
    end
    %%
    
    %Sleeptest trials
    for nT=1:length([trials_st(4).id])
        
        try
            %select N2 and REM
            nSta=find(cellfun(@(x) ismember(trials_st(4).tagSleepScoring(nT),x),Numstages(2:end)))+1;
            if isempty(nSta)
                continue
            end
            
            count(nS,nSta)=count(nS,nSta)+1;
            fprintf('... newSeg in %s subject %s \n',NameStages{nSta}, SubID);
            
            window_time=trials_st(4).trials(nT,1):trials_st(4).trials(nT,2);
            
            for bandf=1:length(bandfreq)
                freqrange1=bandfreq(bandf,:);
                freq_stats{nS,nSta}(count(nS,nSta),bandf,chanSelec) = bandpower(data_full(chanSelec,window_time)',Fs,freqrange1);
            end
            nE_id{nS,nSta}(count(nS,nSta))=nE;
            if strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'L')
                side_Tale{nS,nSta}(count(nS,nSta))=1;
            elseif strcmp(SubjectBehavData.TrialsCaracs(nE).tSide,'R')
                side_Tale{nS,nSta}(count(nS,nSta))=2;
            end
            
            %30s-60s
            if ~isempty(begS)
                for nSeg=1:length(begS)
                    window_time=[trials_st(4).trials(nT,1)+(nSeg-1)*30*Fs]:min([trials_st(4).trials(nT,1)+(nSeg)*30*Fs,trials_st(4).trials(nT,2)]);
                    for bandf=1:length(bandfreq)
                        freqrange1=bandfreq(bandf,:);
                        freq_stats_Seg{nS,nSta}(count(nS,nSta),nSeg,bandf,chanSelec) = bandpower(data_full(chanSelec,window_time)',Fs,bandfreq(bandf,:));
                    end
                    
                end
            end
        catch
            problem_nE=[problem_nE;nS,nE]
        end
        
    end
end
a=datetime;
save([savePath,filesep,userID,'freq_stats',datestr(a) '.mat'],'freq_stats','freq_stats_Seg','begS','endS','subject_id','count','Numstages','NameStages','side_Tale','nE_id','problem_nE','bandfreq')%% plot for subjects

%%% plot Results

%%
goodsubj=cell(1,length(Numstages));
Rmat_freq=[];Rmat_freq_Seg=[];
freq_stats_mean=zeros(length(subject_id),length(Numstages),length(bandfreq),length(chanSelec));
freq_stats_Seg_mean=zeros(length(subject_id),length(Numstages),length(begS),length(bandfreq),length(chanSelec));

for nSta=1:length(Numstages)
    for nS=1:length(subject_id)
        if count(nS,nSta)>0
            for bandf=1:length(bandfreq)
                for channelS=1:length(chanSelec)
                    freq_stats_mean(nS,nSta,bandf,channelS)=mean(freq_stats{nS,nSta}(:,bandf,channelS),1);
                    
                    %                 vecone=ones(count(nS,nSta),1);
                    %                 freq_trial=[freq_stats{nS,nSta}(:,bandf,channelS),nSta*vecone,bandf*vecone,chanSelec(channelS)*vecone,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
                    %                 Rmat_freq=[Rmat_freq;freq_trial];
                    
                    if ~isempty(begS)
                        for nSeg=1:length(begS)
                            freq_stats_Seg_mean(nS,nSta,nSeg,bandf,channelS)=mean(freq_stats_Seg{nS,nSta}(:,nSeg,bandf,channelS),1);
                            %                         freq_trial=[freq_stats_Seg{nS,nSta}(:,nSeg,bandf,channelS),nSeg*vecone,nSta*vecone,bandf*vecone,chanSelec(channelS)*vecone,nS*vecone,nE_id{nS,nSta},side_Tale{nS,nSta}];
                        end
                    end
                end
            end
        end
    end
    goodsubj{nSta}=find(count(:,nSta)>4);
end
%% plot figure1
channelSelect=[65];
plotFigureFreq
