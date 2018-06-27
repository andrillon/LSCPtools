function fname_spm=core_detectSleepRhythms(param)

% load D
D=spm_eeg_load(param.fname_spm);
fname=D.fname; temp=findstr(fname,'_'); 
if ~isempty(temp)
fname(1:temp(end))=[]; 
end
if fname(1)=='M'
    fname(1)=[];
end
fname(findstr(fname,'.mat'):end)=[];

% Check scoring
if exist([param.sleepScoring.Dir filesep param.sleepScoring.Prefix fname '.mat'])~=0
    SleepScoring_filename=[param.sleepScoring.Dir filesep param.sleepScoring.Prefix fname '.mat'];
else
    warning(sprintf('... no scoring for %s!\n',fname))
    SleepScoring_filename=[];
end

% Detecting Slow Waves
paramDetection=param.SWdetection;
[allWaves, slowWaves] = SWsDetectionAlgorithm_forSPM(D, SleepScoring_filename, paramDetection);
% save results
saveName=[param.spm_datapath filesep 'SWdetection_' fname];
save(saveName,'allWaves','slowWaves');

% Detecting Sleep Spindles
paramDetection=param.SpindlesDetection;
[spindles spindlesOutEnergy spindlesOutDuration] = SpindlesDetectionAlgorithm_forSPM(D, SleepScoring_filename, paramDetection);
saveName=[param.spm_datapath filesep 'SPdetection_' fname];
save(saveName,'spindles','spindlesOutEnergy','spindlesOutDuration');

fname_spm=param.fname_spm;
