%%% --------------------------------------------
%%%
%%% Sleep Spindles detection algorithm
%%%
%%% From Andrillon, Nir at al. J Neurosc. 2011
%%%
%%% Description:
%%% Simplified version of the spindle detection algorithm
%%% described in Andrillon, Nir et al. J Neurosc 2011. T
%%%
%%% Input:
%%%     - EEGdata: EEG data, vector
%%%     - SR: sampling rate of the data
%%%     - SleepStages: cell structure with (1) the size of the window used
%%%     for scoring (2) the scoring of each epoch. If empty, threshold will
%%%     be computed overall
%%%     - paramDetection: structure setting the detection parameters.
%%%     Can be left empty to run on default parameters
%%%         -- paramDetection.spindleBand: specify the frequency band used
%%%         to bandpass the EEG and the filter order (default: [11 16 4])
%%%         -- paramDetection.minMax_SD_threshold: 3-scalars vector
%%%         setting the detection thresholds in SD. The first scalar is
%%%         the detection threshold (default: 3), the second the rejection
%%%         threshold (default: 10) and the third is used to define begin/
%%%         end of spindles (default: 1);
%%%         -- paramDetection.minMaxDurations_ms: 2-scalars vector
%%%         setting the minimum (default: 500ms) and maximum (default:
%%%         3000ms) durations of spindles in ms;
%%%         -- paramDetection.minFreqValue: minimum frequency value in Hz
%%%         used to compute spindle frequency (default: 9);
%%%         -- paramDetection.maxFreqValue: maximum frequency value in Hz
%%%         used to compute spindle frequency (default: 17);
%%%         -- paramDetection.PlotFlag: to plot detected spindles or not
%%%         (default: 0, no plot);
%%%
%%% Output:
%%%     - spindles: matrix bearing information on detected spindles:
%%% [SpindleStartTime SpindleEndTime AmplitudePeakTime AmplitudePeakEnergy SpindleFrequency SpindleDuration SleepStage];
%%%
%%%
%%% 21 June 2013 - contact: thomas.andrillon@ens.fr
%%% --------------------------------------------


function [spindles spindlesOutEnergy spindlesOutDuration] = SpindlesDetectionAlgorithm_fieldtrip(data,hdr, SleepScoring_vector, paramDetection)

spindles=[];
spindlesOutEnergy=[];
spindlesOutDuration=[];

SR=hdr.Fs;

if nargin<3
    paramDetection=[];
    fprintf('\n *** Parameters set to default ***\n')
end
if isfield(paramDetection,'spindleBand'),         spindleBand=paramDetection.spindleBand;                 else spindleBand=[11 16 4]; end
if isfield(paramDetection,'minMax_SD_threshold'), minMax_SD_threshold=paramDetection.minMax_SD_threshold; else minMax_SD_threshold=[3 10 1]; end % 3 10 1
if isfield(paramDetection,'minMaxDurations_ms'),  minMaxDurations_ms=paramDetection.minMaxDurations_ms;   else minMaxDurations_ms=[500 2500]; end
if isfield(paramDetection,'minFreqValue'),        minFreqValue=paramDetection.minFreqValue;               else minFreqValue=11; end
if isfield(paramDetection,'maxFreqValue'),        maxFreqValue=paramDetection.maxFreqValue;               else maxFreqValue=17; end
if isfield(paramDetection,'PlotFlag'),            PlotFlag=paramDetection.PlotFlag;                       else PlotFlag=0; end
if isfield(paramDetection,'ChannelSelection'),    ChannelSelection=paramDetection.ChannelSelection;       else ChannelSelection=(1:hdr.nChans); end
if isfield(paramDetection,'NREMonlyForThr'),      NREMonlyForThr=paramDetection.NREMonlyForThr;           else NREMonlyForThr=1; end

% S=[];
% S.D=D;
% S.filter.type='but';
% S.filter.order=spindleBand(3);
% S.filter.band='bandpass';
% S.filter.PHz=spindleBand(1:2);
% S.filter.dir='twopass';
% S.save=0;
% Dfilt=spm_eeg_filter(S);

countChan=0;
for nChan=ChannelSelection
    countChan=countChan+1;
    fprintf('... detecting spindles on %s (%g/%g)',hdr.label{nChan},countChan,length(ChannelSelection))
    if ndims(data)==2
        EEGdata=squeeze(data(nChan,:));
    else
        EEGdata=squeeze(data(nChan,:,:));
        EEGdata=reshape(EEGdata',1,numel(EEGdata));
    end
    %% Filter the EEG within the spindle range
    EEGdata_BP_spindles=bandpass(EEGdata, SR, spindleBand(1), spindleBand(2), spindleBand(3));
    EEGdata_BP_broad=bandpass(EEGdata, SR, 0.1, 30, spindleBand(3));
    EEGdata_BP_high=bandpass(EEGdata, SR, 30, 50, spindleBand(3));
    EEGdata_BP_alpha=bandpass(EEGdata, SR, 8, 10, spindleBand(3));
    %% Select only sleep segments
    if isempty(SleepScoring_vector)
        EpochsOfInterest=logical(ones(1,length(EEGdata)));
        SleepStages_ts=[];
    else
        if length(SleepScoring_vector)~=size(data,2)
            error('Wrong size scoring');
        end
        SleepStages_ts=repmat(SleepScoring_vector',1,size(data,3));
        SleepStages_ts=reshape(SleepStages_ts',1,numel(SleepStages_ts));
        if NREMonlyForThr
            NREM23epochs = find( ismember(SleepStages_ts,[1 2 3])); % Here stages are coded as follows: Wake: 0, N1: 1, N2: 2, N3: 3, REM: 5
            EpochsOfInterest=zeros(1,length(EEGdata));
            EpochsOfInterest(NREM23epochs)=1;
        else
            EpochsOfInterest=logical(ones(1,length(EEGdata)));
        end
    end
    artefactsEpochs=abs(hilbert(EEGdata_BP_high))>30 | abs(EEGdata)>500;
    beg_artefactsEpochs=find(diff(artefactsEpochs)==1);
    end_artefactsEpochs=find(diff(artefactsEpochs)==-1);
    if length(beg_artefactsEpochs)==length(end_artefactsEpochs)+1
        beg_artefactsEpochs(end)=[];
    elseif length(beg_artefactsEpochs)==length(end_artefactsEpochs)-1
        end_artefactsEpochs(1)=[];
    end
    for nE=1:length(beg_artefactsEpochs)
        EpochsOfInterest(max(beg_artefactsEpochs(nE)-SR,1):min(end_artefactsEpochs(nE)+SR,length(EEGdata)))=false;
    end
    
    %% Compute envelope and thresholds
    envelope = (abs(hilbert(EEGdata_BP_spindles)));
    
    envelope_EpochsOfInterest = envelope;
    envelope_EpochsOfInterest(EpochsOfInterest==0)=NaN;
    
    detectionThresholds(1)=nanmedian(envelope_EpochsOfInterest);
    detectionThresholds(2)=nanstd(envelope_EpochsOfInterest);
    
    DetectionThreshold  = detectionThresholds(1) + detectionThresholds(2)*minMax_SD_threshold(1);
    RejectThreshold     = detectionThresholds(1) + detectionThresholds(2)*minMax_SD_threshold(2);
    StartEndThreshold   = detectionThresholds(1) + detectionThresholds(2)*minMax_SD_threshold(3);
    fprintf(' ... detection threshold %g | start/end threshold  %g | rejection threshold %g',DetectionThreshold,StartEndThreshold,RejectThreshold)
    
    %% Initialization of variables to save
    these_spindles = [];
    these_spindlesOutDuration = [];
    these_spindlesOutEnergy = [];
    
    if isnan(detectionThresholds(1))
        return;
    end
    
    %% %% Detection
    % Transform envelope so that every point above the threshold is a 'zero crossing'
    envelopeMinusThreshold = envelope - DetectionThreshold;
    pos_index=zeros(length(envelopeMinusThreshold),1);
    pos_index(find(envelopeMinusThreshold>0))=1; %index of all positive points for EEG
    difference=diff(pos_index); poscross=find(difference==1) ; negcross=find(difference==-1); %find neg ZX and pos ZX
    peaks = find(imregionalmax(envelopeMinusThreshold) == 1);
    peaks(envelopeMinusThreshold(peaks)<0)=[];
    
    if isempty(negcross) || isempty(poscross)
        spindles{nChan}=[];
        spindlesOutDuration{nChan}=[];
        spindlesOutEnergy{nChan}=[];
        fprintf('... no spindle found\n')
        continue;
    end
    
    %% Spindles Between 2 segments
    % Correct anomalies in negcross and poscross
    if max(isempty(poscross),isempty(negcross))==0
        if poscross(end)>negcross(end)
            poscross(end)=[];
        elseif negcross(1)<poscross(1)
            negcross(1)=[];
        end
    end
    if max(isempty(poscross),isempty(negcross))==0
        if (length(poscross) ~= length(negcross)) % this can happen if the spindle occurs when we transition from one segment to the next..
            %         display('WARNING : Spindle between 2 segments')
            if negcross(1)<poscross(1)
                negcross(1)=[];
            elseif poscross(end)>negcross(end)
                poscross(end)=[];
            end
        end
    else
        continue;
    end
    % supress poscross too close from edge (2s)
    poscross(poscross<2*SR)=[];
    poscross(poscross>length(EEGdata)-2*SR)=[];
    
    % Transform envelope relative to StartEnd threshold to compute spindle durations
    envelopeMinusBottomThreshold = envelope - StartEndThreshold;
    pos_index=zeros(length(envelopeMinusBottomThreshold),1);
    pos_index(find(envelopeMinusBottomThreshold>0))=1; %index of all positive points for EEG
    difference=diff(pos_index); poscrossMinimal=find(difference==1) ; negcrossMinimal=find(difference==-1); %find neg ZX and pos ZX
    
    %% Compute spindles candidates properties
    spindleStartTimeArray = []; % we build these arrays as a tool to avoid detecting the same spindle twice
    spindleEndTimeArray   = [];
    sndx=0;
    while sndx<length(poscross)
        
        sndx=sndx+1;
        if isempty(poscrossMinimal) || isempty(negcrossMinimal) % Signal crossed the detection threshold but did not corss twice the StartEnd thr
            continue;
        end
      
        % Compute Spindle Start and End time
        spindleStartTime = findNearestPointBefore(poscross(sndx), poscrossMinimal);
        spindleEndTime   = findNearestPointAfter(negcross(sndx),  negcrossMinimal) ;
        if isempty(spindleStartTime) || isempty(spindleEndTime) || spindleStartTime-SR<1 || spindleEndTime+SR>length(EEGdata_BP_spindles)
            continue;
        end
        if (length(spindleStartTime) .* length(spindleEndTime) == 0) % failed to get start time or end time for some reason
            continue;
        end
        spindleDuration  = spindleEndTime - spindleStartTime + 1;
        
        % Check duration is over 0.5s
        if spindleDuration/SR<0.5
            continue;
        end
        
        % Check that this spindle candidate is not too close from the next one
        % (<500ms). If it is the case, spindles candidates are concatenated
        if sndx~=length(poscross)
            partDuration=spindleDuration;
            while (sndx+1<length(poscross)+1 && isempty(findNearestPointAfter(negcross(sndx+1),  negcrossMinimal))==0 && (findNearestPointBefore(poscross(sndx+1), poscrossMinimal)-spindleEndTime)<0.3*SR)
                spindleStartTime = findNearestPointBefore(poscross(sndx), poscrossMinimal);
                spindleEndTime   = findNearestPointAfter(negcross(sndx+1),  negcrossMinimal) ;
                AddpartDuration=findNearestPointAfter(negcross(sndx+1),  negcrossMinimal)-findNearestPointBefore(poscross(sndx+1), poscrossMinimal);
                partDuration  = partDuration + AddpartDuration;
                poscross(sndx+1)=[];
                negcross(sndx)=[];
            end
        end
         if spindleStartTime-SR<1 || spindleEndTime+SR>length(EEGdata_BP_spindles)
            continue;
         end
        
        % Compute Spindle Duration
        spindleDuration  = spindleEndTime - spindleStartTime + 1;
        
        % Make sure this is not just another peak for the same spindle (as judged by its start and end times) if so continue and ignore..
        if (size(intersect([spindleStartTime], spindleStartTimeArray) > 0))
            continue;
        elseif (size(intersect([spindleEndTime], spindleEndTimeArray) > 0))
            continue;
        end
        spindleStartTimeArray = [spindleStartTimeArray; spindleStartTime];
        spindleEndTimeArray = [spindleEndTimeArray; spindleEndTime];
        
        % Compute Peak in amplitude time and power
        relevantPeakIdx = find ( (peaks > spindleStartTime) & (peaks < spindleEndTime));
        peakTimes = peaks(relevantPeakIdx);
        peakWithMaxValue = find(envelope(peakTimes) == max(envelope(peakTimes)));
        peakTime = peakTimes(peakWithMaxValue);
        if isempty(peakTime), continue; end;
        peakTime=peakTime(1);
        peakEnergy = envelope(peakTime);
        peakEnergyNorm = (envelope(peakTime)-detectionThresholds(1))/detectionThresholds(2);
        
        % discard spindle if spindles occurs during artefacted epochs
        if max(abs(hilbert(EEGdata_BP_high(spindleStartTime-SR:spindleEndTime+SR))))>30
            continue;
        end
        %         % Power Spectrum
        %         [faxis,myPSDs]= get_PowerSpec_new(EEGdata(max(peakTime-499,1):min(peakTime+500,10000)),1/1000,1);
        %         faxis=interp1((1:length(faxis)),faxis,(1:1/10:length(faxis)),'linear');
        %         myPSDs=interp1((1:length(myPSDs)),myPSDs,(1:1/10:length(myPSDs)),'linear');
        
        % Compute Spindle candidate frequency using the power spectrum
        % (precision: 1Hz)
        [faxisFilt,myPSDsFilt]=get_PowerSpec_new(EEGdata_BP_spindles(peakTime-0.5*SR:peakTime+0.5*SR),1/SR,length(-0.5*SR:0.5*SR)/SR);
        faxisFilt=interp1((1:length(faxisFilt)),faxisFilt,(1:1/10:length(faxisFilt)),'linear');
        myPSDsFilt=interp1((1:length(myPSDsFilt)),myPSDsFilt,(1:1/10:length(myPSDsFilt)),'linear');
        freqSpindle=faxisFilt(find(myPSDsFilt==max(myPSDsFilt(find(faxisFilt>=minFreqValue & faxisFilt<=maxFreqValue)))));
        
        %%% compute ratio with alpha power
        PowerSP=mean(abs(EEGdata_BP_spindles(spindleStartTime:spindleEndTime)));
        PowerAlpha=mean(abs(EEGdata_BP_alpha(spindleStartTime:spindleEndTime)));
        
        % Stage
        if ~isempty(SleepStages_ts)
            currentStage=min(SleepStages_ts(spindleStartTime:spindleEndTime));
        else
            currentStage=NaN;
        end
        % Sum-up
        currentSpindle=[spindleStartTime spindleEndTime peakTime peakEnergy peakEnergyNorm freqSpindle spindleDuration PowerSP PowerAlpha currentStage];
        
        %% %% Sort Spindles
        % Candidate discarded: Peak Energy too highg
        if (peakEnergy > RejectThreshold)
            these_spindlesOutEnergy = [these_spindlesOutEnergy ; currentSpindle];
            continue;
        end
        
        % Candidate discarded: Out Because of  Duration
        if ( (spindleDuration > minMaxDurations_ms(1)*SR/1000) && (spindleDuration < minMaxDurations_ms(2)*SR/1000) )
            these_spindles = [these_spindles ; currentSpindle];
        else
            these_spindlesOutDuration = [these_spindlesOutDuration ; currentSpindle];
            continue;
        end
        
        
    end
    spindles{nChan}=these_spindles;
    spindlesOutDuration{nChan}=these_spindlesOutDuration;
    spindlesOutEnergy{nChan}=these_spindlesOutEnergy;
    
    fprintf('... %g spindles detected',size(spindles{nChan},1))
    fprintf('\n')
    %% %% Plotting
    if (PlotFlag==1) % plot detected spindles
        for nSpin=1:size(these_spindles,1)
            % Set figure parameters
            f1 = figure;
            xAxis = -5:1/SR:5;
            xPoints=-5*SR:5*SR;
            set(f1, 'Position',get(0,'ScreenSize'));
            
            subplot(4,1,1);
            plot(xAxis, (EEGdata_BP_broad(these_spindles(nSpin,3)+xPoints)), 'b'); hold on; % axis([1 10 -10 10]);
            
            subplot(4,1,2);
            plot(xAxis, (EEGdata_BP_spindles(these_spindles(nSpin,3)+xPoints)), 'r:'); hold on;
            plot(xAxis, envelope(these_spindles(nSpin,3)+xPoints), 'k');   hold on; myAxis = axis(); hold on;
            line([myAxis(1,1) myAxis(1,2)], [DetectionThreshold DetectionThreshold], 'Color', 'k', 'LineWidth', 1);hold on;
            line([myAxis(1,1) myAxis(1,2)], [RejectThreshold RejectThreshold], 'Color', 'b', 'LineWidth', 1);  hold on;
            line([myAxis(1,1) myAxis(1,2)], [StartEndThreshold StartEndThreshold], 'Color', 'g', 'LineWidth', 1); hold on;
            hold on;
            plot(0, envelope(these_spindles(nSpin,3)), 'or', 'MarkerFaceColor',[1 0 0], 'MarkerSize',7); hold on;
            plot((these_spindles(nSpin,1)-these_spindles(nSpin,3))/SR, envelope(these_spindles(nSpin,1)), 'oc', 'MarkerFaceColor',[0 1 1], 'MarkerSize',7); hold on;
            plot((these_spindles(nSpin,2)-these_spindles(nSpin,3))/SR, envelope(these_spindles(nSpin,2)), 'og', 'MarkerFaceColor',[0 1 0], 'MarkerSize',7); hold on;
            
            axis([min(xAxis) max(xAxis) myAxis(1,3) RejectThreshold*1.2]);
            title(sprintf('Duration: %g, Frequency: %g, Stage: %g', these_spindles(nSpin,6)/SR, these_spindles(nSpin,5), these_spindles(nSpin,7)));
            
            subplot(4,1,3);
            plot(xAxis, (EEGdata_BP_high(these_spindles(nSpin,3)+xPoints)), 'b'); hold on; % axis([1 10 -10 10]);
            
            subplot(4,1,4)
            [faxisFilt, myPSDsFilt]=get_PowerSpec_new(EEGdata(these_spindles(nSpin,2):these_spindles(nSpin,3)),1/SR,(these_spindles(nSpin,3)-these_spindles(nSpin,2))/SR,1,0);
            plot(faxisFilt,myPSDsFilt); xlim([0 30])
            %         line([minFreqValue minFreqValue], [min(myPSDsFilt) max(myPSDsFilt)], 'Color','r')
            %         line([maxFreqValue maxFreqValue], [min(myPSDsFilt) max(myPSDsFilt)], 'Color','r')
            axis('tight'); xlim([0 30])
            pause;
            close(f1);
        end
    end
end
%
%
%
%
% %% %% Spindles quick analysis
% % 11 : W / 12 : REM / 13 : N1 / 14 : N2 / 15 : N3
% if isempty(these_spindles)==0
%
%     numSpindlesIn_NREM2 = length(find( these_spindles(:, end) == 14));
%     numSpindlesIn_NREM3 = length(find( these_spindles(:, end) == 15));
%     numSpindlesIn_NREM = length(find( these_spindles(:, end) >= 14));
%     numSpindlesIn_REM  = length(find( these_spindles(:, end) == 12));
%
%     length_NREM2_min  = length(find(SleepStages_ts == 14)) / 6;
%     length_NREM3_min  = length(find(SleepStages_ts == 15)) / 6;
%     length_NREM_min = length(find(SleepStages_ts >= 14)) / 6;
%     length_REM_min  = length(find(SleepStages_ts == 12)) / 6;
%
%     spindlesPerMin_NREM2 = numSpindlesIn_NREM2 / length_NREM2_min;
%     spindlesPerMin_NREM3 = numSpindlesIn_NREM3 / length_NREM3_min;
%     spindlesPerMin_NREM = numSpindlesIn_NREM / length_NREM_min;
%     spindlesPerMin_REM  = numSpindlesIn_REM  / length_REM_min;
%     % fprintf('Working on %s, found %d spindles, %d percent in NREM, %d percent in REM\n', EEG_filename_spindles, size(these_spindles,1), percentageOfSpindlesIn_NREM, percentageOfSpindlesIn_REM);
%     fprintf('Found %d spindles, %.2f spindles/min in NREM, %.2f spindles/min in REM\n', size(these_spindles,1), spindlesPerMin_NREM, spindlesPerMin_REM);
%     fprintf('%.2f spindles/min in NREM2, %.2f spindles/min in NREM3\n', spindlesPerMin_NREM2, spindlesPerMin_NREM3);
%
% end


%% %% Save Files
% save what you want to keep

end

%%%%%%%%%% Embedded functions
function BP = bandpass(timecourse, SamplingRate, low_cut, high_cut, filterOrder)
if (nargin < 5)
    filterOrder = 2;
end
[b, a] = butter(filterOrder, [(low_cut/SamplingRate)*2 (high_cut/SamplingRate)*2]);
BP = filtfilt(b, a, timecourse );
end

function  nearestPoint = findNearestPointBefore(seedTime, array)
nearestPoint = [];
timeDiffAll = seedTime - array;
idxOfTimeDiffBefore = find(timeDiffAll > 0);
positiveTimeDiffBefore = timeDiffAll(idxOfTimeDiffBefore);
idxOfNearestPoint = find(positiveTimeDiffBefore == min(positiveTimeDiffBefore));
minimalDiffereceBefore = positiveTimeDiffBefore(idxOfNearestPoint);
nearestPoint = seedTime - minimalDiffereceBefore;
end

function  nearestPoint = findNearestPointAfter(seedTime, array)
nearestPoint = [];
timeDiffAll = array - seedTime;
idxOfTimeDiffAfter = find(timeDiffAll > 0);
positiveTimeDiffAfter = timeDiffAll(idxOfTimeDiffAfter);
idxOfNearestPoint = find(positiveTimeDiffAfter == min(positiveTimeDiffAfter));
minimalDiffereceAfter = positiveTimeDiffAfter(idxOfNearestPoint);
nearestPoint = seedTime + minimalDiffereceAfter;
end

function [faxis, pow] = get_PowerSpec_new(signal, SamplingRate, lengthSec, DecibelsFlag ,plotFlag)
% For example for a 10sec segment sampled at millisec resolution use:
% [faxis, pow] = get_PowerSpec_new(mySignal, 1/1000, 10);
if (nargin < 5), plotFlag = 0; end
if (nargin < 4), DecibelsFlag = 0; end
% get power
pow = (abs(fft(signal)).^2)/length(signal);
% convert to decibels
if DecibelsFlag==1
    pow = 10*log10(pow/max(pow));
end
% first half of data without negative frequencies
pow = pow(1:min(floor(length(signal)/2)+1,length(signal)));
% define df and fNQ
df = 1/lengthSec;
fNQ = 1/SamplingRate/2;

faxis = (0:df:fNQ);
if (plotFlag)
    plot(faxis, pow);
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
end
end
