% Function to compute LRPs

function LRD=compute_LRD3(D,param)

% DEBUG
if ~isfield(param,'chanLabels')
    chanLabels.left={'C3-M'};
    chanLabels.right={'C4-M'};
else
    chanLabels=param.chanLabels;
end
if ~isfield(param,'respLabels')
    respLabels.left='LeftCorrect';
    respLabels.right='RightCorrect';
else
    respLabels=param.respLabels;
end
if ~isfield(param,'badT'), badT=logical(zeros(1,D.ntrials)); else badT=param.badT; end
if ~isfield(param,'minT'), minT=0; else minT=param.minT; end
if ~isfield(param,'avE'), avE=1; else avE=param.avE; end
if ~isfield(param,'freqB'), freqB=D.frequencies; else freqB=param.freqB; end
if ~isfield(param,'dbFlag'), dbFlag=0; else dbFlag=param.dbFlag; end

C3_idx=[];
if iscell(chanLabels.left)
    for nC=1:length(chanLabels.left)
        C3_idx=[C3_idx match_str(D.chanlabels,chanLabels.left{nC})];
    end
else
    C3_idx=find(cellfun(@any, regexp(D.conditions, chanLabels.left)));
end
C4_idx=[];
if iscell(chanLabels.right)
    for nC=1:length(chanLabels.right)
        C4_idx=[C4_idx match_str(D.chanlabels,chanLabels.right{nC})];
    end
else
    C4_idx=find(cellfun(@any, regexp(D.conditions, chanLabels.right)));
end

leftResp_idx=find(cellfun(@any, regexp(D.conditions,respLabels.left)));
rightResp_idx=find(cellfun(@any, regexp(D.conditions,respLabels.right)));

% Clean data
orileftResp_idx=leftResp_idx;
orirightResp_idx=rightResp_idx;
if ~isempty(badT)
    myN=D.fname;
    D2=spm_eeg_load([D.path filesep myN(4:end)]);
    % raw cleaning
    badT_abs=(max(abs(squeeze(mean(D2(C3_idx,:,:),1))))>badT(1) | max(abs(squeeze(mean(D2(C4_idx,:,:),1))))>badT(1));
%     badT_gra=(max(abs(squeeze(mean(D2(C3_idx,1:size(D,2)-39,:),1))-squeeze(mean(D2(C3_idx,40:size(D2,2),:),1))))>badT(2) | max(abs(squeeze(mean(D2(C4_idx,1:size(D2,2)-39,:),1))-squeeze(mean(D2(C4_idx,40:size(D2,2),:),1))))>badT(2));
    badTr=(badT_abs); % | badT_gra);
    
    % lEOG=squeeze(data(match_str(D.chanlabels,'lEOG'),:,:));
    % lEOG=(lEOG-mean(reshape(lEOG,1,numel(lEOG))))/std(reshape(lEOG,1,numel(lEOG)));
    % rEOG=squeeze(data(match_str(D.chanlabels,'rEOG'),:,:));
    % rEOG=(rEOG-mean(reshape(rEOG,1,numel(rEOG))))/std(reshape(rEOG,1,numel(rEOG)));
    % badT_eog=(max(abs(lEOG),[],1)>badT(1) | max(abs(rEOG),[],1)>badT(1));
    %
    % chinEMG=squeeze(data(match_str(D.chanlabels,'chinEMG'),:,:));
    % chinEMG=(chinEMG-mean(reshape(chinEMG,1,numel(chinEMG))))/std(reshape(chinEMG,1,numel(chinEMG)));
    % badT_emg=max(abs(chinEMG),[],1)>badT(1);
    % badTr=(badT_eog | badT_emg);
    
    [c toRl c2]=intersect(leftResp_idx,find(badTr));
    [c toRr c2]=intersect(rightResp_idx,find(badTr));
    
    leftResp_idx(toRl)=[];
    rightResp_idx(toRr)=[];
end
if length(leftResp_idx)<minT ||length(rightResp_idx)<minT
    fprintf('... right cond: %g (%g) | left cond: %g (%g)\n',length(leftResp_idx),length(leftResp_idx)/length(orileftResp_idx)*100,length(rightResp_idx),length(rightResp_idx)/length(orirightResp_idx)*100)
    fprintf('... not enough trials for at least one condition (<%g)\n',minT)
    LRD = nan(1,size(D,2),size(D,3));
else
    fprintf('... right cond: %g (%g) | left cond: %g (%g)\n',length(leftResp_idx),length(leftResp_idx)/length(orileftResp_idx)*100,length(rightResp_idx),length(rightResp_idx)/length(orirightResp_idx)*100)
    
    if dbFlag
        data=D(:,:,:,:);
        temp1=repmat((mean(data(:,:,D.time<0,:),3)),[1 1 size(D,3) 1]);
        temp2=(data-temp1)./temp1*100;
        data=(temp2);
    else
                data=D(:,:,:,:);
    end
    % LRP formula
    % LRP= 0.5 [meanleft hand response(C4V–C3V) + meanright hand response(C3V–C4V)]
            [usl1 myF usl2]=intersect(D.frequencies,freqB);
    if avE
        LRD = squeeze(0.5 * (...
            nanmean(nanmean(data(C4_idx,myF,:,leftResp_idx)-data(C3_idx,myF,:,leftResp_idx),1),4) + ...
            nanmean(nanmean(data(C3_idx,myF,:,rightResp_idx)-data(C4_idx,myF,:,rightResp_idx),1),4)));
    else
        LRD = squeeze(0.5 * (...
            nanmean(data(C4_idx,myF,:,leftResp_idx)-data(C3_idx,myF,:,leftResp_idx),4) + ...
            nanmean(data(C3_idx,myF,:,rightResp_idx)-data(C4_idx,myF,:,rightResp_idx),4)));
    end
    
    
end