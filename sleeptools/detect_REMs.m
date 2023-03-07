%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function detect_REMs.m
%%% The aim of this function is to detect Rapid-Eye Movents in sleep EOG data
%%% The function is based on the asumptions that REMs appears on EOG1 and EOG2 as big anti-phase delfections
%%%
%%% Input:
%%% - eog1, eog2:      vectors of EOG data (1xtime stamps) (!!! already
%%% re-referenced to opposite MASTOID)
%%% - param: 
%%%     - Fs:              sampling rate
%%%     - thresholdParam:  [SD_BP SD_high] SD threshold for bandpass and high-bandpass filtered data
%%%     - exclusionParam:  exckusion parameters: (1) max time duration
%%%     between negcross and poscross (2) minimal duration between REMs (in
%%%     seconds), (3) minum slope
%%%     - displayFlag:     set to 1 if you want to display summary
%%%
%%% Output:
%%% - REMs:            Detected REMs. The matrix stores the following informations:
%%%                   (1) Segment (2) Begin (3) End (4) Peak time (5) slope peak time (6) Peak EOG times (7) Peak ampl EOG1 (8) Peak ampl EOG2 (9) Slope max 1 (10) Slope max 2 (11) Duration (12) Stage
%%% - false_detection: False detected events
%%% - thresholds     : thesholds used for (1) product and (2) sum
%%%
%%% Author: Thomas Andrillon (thomas.andrillon@ens.fr or thomas.andrillon@gmail.com)
%%% Created: 2-14-12
%%%
%%% History:
%%%     v1: Initial detection with a series of criterion
%%%     v2: Hopefully cleaned detection
%%% Last update: 2-28-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [REMs , false_detection] = detect_REMs(eog1, eog2, param)

%% Initialization
REMs=[];
false_detection=[];
if ~isfield(param,'displayFlag')
    displayFlag=0;
else
    displayFlag=param.displayFlag;
end
thresholdParam=param.thresholdParam;
exclusionParam=param.exclusionParam;

%% Load and reshape data
% fprintf('>>> preparing data\n')
% tic;
Fs=param.Fs;

% (1) Broad bandpass
eog1_Broad=bandpass(eog1,Fs,0.1,30,3);
eog2_Broad=bandpass(eog2,Fs,0.1,30,3);

% (2) BP 3 Hz
eog1_BP=bandpass(eog1,Fs,0.1,3,3);
eog2_BP=bandpass(eog2,Fs,0.1,3,3);

% (3) Zscore
eog1_BP   = zscore(eog1_BP);
eog2_BP   = zscore(eog2_BP);

% (4) gradient
grad_eog1=abs(eog1_BP(:,round(81*Fs/1000):end)-eog1_BP(:,1:end-round(80*Fs/1000)));
grad_eog2=abs(eog2_BP(:,round(81*Fs/1000):end)-eog2_BP(:,1:end-round(80*Fs/1000)));

% fprintf('  ... took ... %g  s\n',toc)

%% Detection REMs based on different criterions
% fprintf('\n>>>> Detect REMs on eog1 ...\n')
eog_BP=eog1_BP;
eog_Broad=eog1_Broad;
grad_eog=grad_eog1;

FlagEM=1;
eog_BP_aboveThr=abs(eog_BP)>thresholdParam(1);

% +1 when eog crosses the threshold
pos_candidates=find(diff(eog_BP_aboveThr)==1);
lastREM=-1;

for nREM=1:length(pos_candidates)
    
    % (1) consider one 'positive' cross
    pos_cross=pos_candidates(nREM);
    begin=pos_cross;
    
    % (2) find the next negative cross of threshold if any (-1) before 2s  e.g. or before end of segment
    neg_candidates=pos_cross+find(diff(eog_BP_aboveThr(pos_cross:min(pos_cross+exclusionParam(1)*Fs,length(eog_BP_aboveThr))))==-1)-1;
    if ~isempty(neg_candidates)
        neg_cross=neg_candidates(1);
        ending=neg_cross;
    else
        false_detection=[false_detection ; [begin, ... (2) Begin
            0, ... (3) End
            NaN,...
            FlagEM,...
            1]]; %
        continue;
    end
    
    % (3) find peak between these two crossings
    peak_candidates=find(diff(diff(abs(eog_BP(pos_cross:neg_cross)))>0));
    if length(peak_candidates)==0 % length(peak_candidates)~=1 || length(peak_candidates2)~=1
        false_detection=[false_detection ; [
            begin, ... (2) Begin
            ending, ... (3) End
            NaN,...
            FlagEM,...
            2]]; %
        continue;
    else
        peak_cross=pos_cross+round(peak_candidates(1))-1;
    end
    
    % (4) if last REM peak too close: discard, else update
    if abs(peak_cross-lastREM)<exclusionParam(2)*Fs/4
        %         lastREM=peak_cross;
        false_detection=[false_detection ; [
            begin, ... (2) Begin
            ending, ... (3) End
            eog_BP(peak_cross),...
            FlagEM,...
            3]]; %
        continue;
    end
    
    % (6) find real peak in broadly band-passed EOG: if outside
    % pos_cross:neg_cross, discard
    % get broad BP eog
    if peak_cross-0.5*Fs<1 || peak_cross+0.5*Fs>length(eog1)
        continue;
    end
    thisrealeog=eog_Broad((-0.5*Fs:0.5*Fs)+peak_cross);
    thisrealeog=abs(thisrealeog-nanmean(thisrealeog(1:200)));
    max_real1=max(thisrealeog);
    timemax_real1=find(thisrealeog==max_real1); timemax_real1=timemax_real1(1);
    peakTimeeog=peak_cross - 0.5*Fs + round(timemax_real1);
    
    
    % (7) Find slop maximum
    thisgrad1=grad_eog((round(-0.3*Fs):0)+peak_cross);
    maxslope1=max(thisgrad1);
    timemaxslope1=find(thisgrad1==maxslope1); timemaxslope1=timemaxslope1(1);
    slopeTimeeog=peak_cross - 0.5*Fs + round(timemaxslope1);
    
    if (maxslope1>exclusionParam(3))
        % (1) Segment (2) Begin (3) End (4) Peak time (5) slope peak time (6) Peak EOG times (7) Peak ampl eog (8) Peak ampl EOG2 (9) Slope max 1 (10) Slope max 2 (11) Duration (12) Stage
        lastREM=peak_cross;
        REMs=[REMs ; [
            begin, ... (1) Begin
            ending, ... (2) End
            peak_cross ... (3) Peak prod time
            slopeTimeeog, ... (4) slope peak time
            peakTimeeog, ... (5) Peak EOGs time
            eog_BP(peak_cross), ... (6) Peak ampl eog
            maxslope1, ... (7) Slope max 1
            (neg_cross-pos_cross), ... (8) Duration
            FlagEM]]; %(9) vertical or horizontal
    else
        lastREM=peak_cross;
        false_detection=[false_detection ; [
            begin, ... (2) Begin
            ending, ... (3) End
            eog_BP(peak_cross),...
            FlagEM,...
            4]];
    end   
end

%%
% fprintf('\n>>>> Detect REMs on VEOG ...\n')
eog_BP=eog2_BP;
eog_Broad=eog2_Broad;
grad_eog=grad_eog2;
FlagEM=2;
eog_BP_aboveThr=abs(eog_BP)>thresholdParam(1);

% +1 when eog crosses the threshold
pos_candidates=find(diff(eog_BP_aboveThr)==1);
lastREM=-1;

for nREM=1:length(pos_candidates)
    
    % (1) consider one 'positive' cross
    pos_cross=pos_candidates(nREM);
    begin=pos_cross;
    
    % (2) find the next negative cross of threshold if any (-1) before 2s  e.g. or before end of segment
    neg_candidates=pos_cross+find(diff(eog_BP_aboveThr(pos_cross:min(pos_cross+exclusionParam(1)*Fs,length(eog_BP_aboveThr))))==-1)-1;
    if ~isempty(neg_candidates)
        neg_cross=neg_candidates(1);
        ending=neg_cross;
    else
        false_detection=[false_detection ; [begin, ... (2) Begin
            0, ... (3) End
            NaN,...
            FlagEM,...
            1]]; %
        continue;
    end
    
    % (3) find peak between these two crossings
    peak_candidates=find(diff(diff(abs(eog_BP(pos_cross:neg_cross)))>0));
    if length(peak_candidates)==0 % length(peak_candidates)~=1 || length(peak_candidates2)~=1
        false_detection=[false_detection ; [
            begin, ... (2) Begin
            ending, ... (3) End
            NaN,...
            FlagEM,...
            2]]; %
        continue;
    else
        peak_cross=pos_cross+round(peak_candidates(1))-1;
    end
    
    % (4) if last REM peak too close: discard, else update
    if abs(peak_cross-lastREM)<exclusionParam(2)*Fs/4
        %         lastREM=peak_cross;
        false_detection=[false_detection ; [
            begin, ... (2) Begin
            ending, ... (3) End
            eog_BP(peak_cross),...
            FlagEM,...
            3]]; %
        continue;
    end
    
    % (6) find real peak in broadly band-passed EOG: if outside
    % pos_cross:neg_cross, discard
    % get broad BP eog
    if peak_cross-0.5*Fs<1 || peak_cross+0.5*Fs>length(VEOG)
        continue;
    end
    thisrealeog=eog_Broad((-0.5*Fs:0.5*Fs)+peak_cross);
    thisrealeog=abs(thisrealeog-nanmean(thisrealeog(1:200)));
    max_real1=max(thisrealeog);
    timemax_real1=find(thisrealeog==max_real1); timemax_real1=timemax_real1(1);
    peakTimeeog=peak_cross - 0.5*Fs + round(timemax_real1);
    
    
    % (7) Find slop maximum
    thisgrad1=grad_eog((round(-0.3*Fs):0)+peak_cross);
    maxslope1=max(thisgrad1);
    timemaxslope1=find(thisgrad1==maxslope1); timemaxslope1=timemaxslope1(1);
    slopeTimeeog=peak_cross - 0.5*Fs + round(timemaxslope1);
    
    if (maxslope1>exclusionParam(3))
        % (1) Segment (2) Begin (3) End (4) Peak time (5) slope peak time (6) Peak EOG times (7) Peak ampl eog (8) Peak ampl EOG2 (9) Slope max 1 (10) Slope max 2 (11) Duration (12) Stage
        lastREM=peak_cross;
        REMs=[REMs ; [
            begin, ... (1) Begin
            ending, ... (2) End
            peak_cross ... (3) Peak prod time
            slopeTimeeog, ... (4) slope peak time
            peakTimeeog, ... (5) Peak EOGs time
            eog_BP(peak_cross), ... (6) Peak ampl eog
            maxslope1, ... (7) Slope max 1
            (neg_cross-pos_cross), ... (8) Duration
            FlagEM]]; %(9) vertical or horizontal
    else
        lastREM=peak_cross;
        false_detection=[false_detection ; [
            begin, ... (2) Begin
            ending, ... (3) End
            eog_BP(peak_cross),...
            FlagEM,...
            4]];
    end   
end

fprintf('\n>>>> %g horizontal (%2.2f /min) and %g vertical (%2.2f /min) REMs detected\n\n',sum(REMs(:,end)==1),sum(REMs(:,end)==1)/(length(eog1)/param.Fs/60),sum(REMs(:,end)==2),sum(REMs(:,end)==2)/(length(eog1)/param.Fs/60))

%% Give an overview of results
if displayFlag
    fprintf('>>> %g REMs detected out of %g\n',size(REMs,1),size(REMs,1)+size(false_detection,1))
    %%
    figure;
    plot(eog1_Broad,'r'); hold on;
    plot(eog2_Broad,'b'); hold on;
    % REMs
    scatter(REMs(REMs(:,end)==1,1),eog1_Broad(REMs(REMs(:,end)==1,1)),'og','filled')
    scatter(REMs(REMs(:,end)==2,1),eog2_Broad(REMs(REMs(:,end)==2,1)),'dg','filled')
    % FAs
    TypeMarker='ds+x';
    for nType=1:4
    scatter(false_detection(false_detection(:,end)==nType & false_detection(:,end-1)==1,1),eog1_Broad(false_detection(false_detection(:,end)==nType & false_detection(:,end-1)==1,1)),'Marker',TypeMarker(nType),'MarkerFaceColor','r')
    scatter(false_detection(false_detection(:,end)==nType & false_detection(:,end-1)==2,1),eog2_Broad(false_detection(false_detection(:,end)==nType & false_detection(:,end-1)==2,1)),'Marker',TypeMarker(nType),'MarkerFaceColor','m')
    end
    figure;
    hist(false_detection(:,end),12)
    
    figure; TA_hREMs=[]; TA_vREMs=[];
    REMs(REMs(:,1)<1*Fs+1,:)=[];
    REMs(REMs(:,1)>size(eog1_BP,2)-2*Fs,:)=[];
    for nR=1:size(REMs,1)
        if REMs(nR,end)==1
            TA_hREMs=[TA_hREMs ; eog1_Broad(REMs(nR,1)+(-1*Fs:2*Fs))*sign(REMs(nR,6))-mean(eog1_Broad(REMs(nR,1)+(-1*Fs:-0.5*Fs)))];
            
        else
            TA_vREMs=[TA_vREMs ; eog2_Broad(REMs(nR,1)+(-1*Fs:2*Fs))*sign(REMs(nR,6))-mean(eog2_Broad(REMs(nR,1)+(-1*Fs:-0.5*Fs)))];
        end
%         TA_REMs(2,nR,:)=eog2(REMs(nR,1)+(-1*Fs:2*Fs))*sign(REMs(nR,7));
    end
    plot(TA_hREMs','r'); hold on;
    plot(mean(TA_hREMs),'Color','k','LineWidth',2)
%     plot(squeeze(mean(TA_REMs(2,:,:),2)),'b');
    
    for nFA=1:11
        TA_FAs=[];
        nC=0;
        myfalse_detection=false_detection(false_detection(:,end)==nFA,:);
        for nR2=1:size(myfalse_detection,1)
            if myfalse_detection(nR2,1)>1*Fs && myfalse_detection(nR2,1)<length(eog)-2*Fs
                if nFA<3
                    nC=nC+1;
                    TA_FAs(1,nC,:)=eog(myfalse_detection(nR2,1)+(-1*Fs:2*Fs));
                    TA_FAs(2,nC,:)=eog2(myfalse_detection(nR2,1)+(-1*Fs:2*Fs));
                else
                    nC=nC+1;
                    TA_FAs(1,nC,:)=eog(myfalse_detection(nR2,1)+(-1*Fs:2*Fs))*sign(myfalse_detection(nR2,3));
                    TA_FAs(2,nC,:)=eog2(myfalse_detection(nR2,1)+(-1*Fs:2*Fs))*sign(myfalse_detection(nR2,3));
                end
            end
        end
        if ~isempty(TA_FAs)
            subplot(3,4,nFA+1)
            plot(squeeze(TA_FAs(1,:,:))','r'); hold on;
            plot(squeeze(TA_FAs(2,:,:))','b');
        end
    end
    
end



