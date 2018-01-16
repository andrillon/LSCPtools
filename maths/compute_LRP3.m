% Function to compute LRPs

function LRP=compute_LRP3(D,param)

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
if ~isfield(param,'maxT'), maxT=[]; else maxT=param.maxT; end
if ~isfield(param,'plotFlag'), plotFlag=0; else plotFlag=param.plotFlag; end
if ~isfield(param,'envFlag'), envFlag=0; else envFlag=param.envFlag; end
if ~isfield(param,'reRef'), reRef=[]; else reRef=param.reRef; end
if ~isfield(param,'bandP'), bandP=[]; else bandP=param.bandP; end
if ~isfield(param,'smth'), smth=[]; else smth=param.smth; end
if ~isfield(param,'avE'), avE=1; else avE=param.avE; end

C3_idx=[];
if iscell(chanLabels.left)
    for nC=1:length(chanLabels.left)
        C3_idx=[C3_idx match_str(D.chanlabels,chanLabels.left{nC})];
    end
else
    C3_idx=find(cellfun(@any, regexp(D.chanlabels, chanLabels.left)));
end
C4_idx=[];
if iscell(chanLabels.right)
    for nC=1:length(chanLabels.right)
        C4_idx=[C4_idx match_str(D.chanlabels,chanLabels.right{nC})];
    end
else
    C4_idx=find(cellfun(@any, regexp(D.chanlabels, chanLabels.right)));
end

leftResp_idx=find(cellfun(@any, regexp(param.trialLabels,respLabels.left)));
rightResp_idx=find(cellfun(@any, regexp(param.trialLabels,respLabels.right)));

% Transform data
data=D(:,:,:);
if ~isempty(reRef)
    fprintf('... re-referencing data!\n')
    if strcmp(reRef,'all')
        newRefChan=D.meegchannels;
        data=data-repmat(mean(data(newRefChan,:,:),1),[size(data,1) 1 1]);
    else
        newRefChan=match_str(D.chanlabels,reRef);
        data=data-repmat(mean(data(newRefChan,:,:),1),[size(data,1) 1 1]);
    end
end


% Clean data
orileftResp_idx=leftResp_idx;
orirightResp_idx=rightResp_idx;
if ~isempty(badT)
    % raw cleaning
    badT_abs=(max(abs(squeeze(mean(data(C3_idx,:,:),1))))>badT(1) | max(abs(squeeze(mean(data(C4_idx,:,:),1))))>badT(1));
    badT_abs_bs=(max(abs(squeeze(mean(data(C3_idx,D.time<0,:),1))))>badT(2) | max(abs(squeeze(mean(data(C4_idx,D.time<0,:),1))))>badT(2));
    badTr=(badT_abs | badT_abs_bs);
    
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
if ~isempty(maxT)
    warning(sprintf('TAKING ONLY THE FIRST %g TRIALS!',maxT))
    leftResp_idx(maxT:end)=[];
    rightResp_idx(maxT:end)=[];
end
if length(leftResp_idx)<minT ||length(rightResp_idx)<minT
    fprintf('... right cond: %g (%g) | left cond: %g (%g)\n',length(leftResp_idx),length(leftResp_idx)/length(orileftResp_idx)*100,length(rightResp_idx),length(rightResp_idx)/length(orirightResp_idx)*100)
    fprintf('... not enough trials for at least one condition (<%g)\n',minT)
    LRP = nan(1,size(D,2));
else
    fprintf('... right cond: %g (%g) | left cond: %g (%g)\n',length(leftResp_idx),length(leftResp_idx)/length(orileftResp_idx)*100,length(rightResp_idx),length(rightResp_idx)/length(orirightResp_idx)*100)
    
    if envFlag
        fprintf('... extracting envelope ...\n')
        for nC=[C4_idx C3_idx]
            for nT=[leftResp_idx'  rightResp_idx']
                data(nC,:,nT)=(abs(hilbert(squeeze(data(nC,:,nT)))));
                data(nC,:,nT)=data(nC,:,nT)-mean(data(nC,find(D.time<0),nT),2);
            end
        end
    end
    if ~isempty(bandP)
        fprintf('... bandpassing ...\n')
        for nC=[C4_idx C3_idx]
            for nT=[leftResp_idx' rightResp_idx']
                data(nC,:,nT)=(abs(hilbert(bandpass(squeeze(data(nC,:,nT)),D.fsample,bandP(1),bandP(2),3))));
                data(nC,:,nT)=data(nC,:,nT)-mean(data(nC,find(D.time<0),nT),2);
            end
        end
    end
    if ~isempty(smth)
        fprintf('... smotthing ...\n')
        for nC=[C4_idx C3_idx]
            for nT=[leftResp_idx rightResp_idx]
                data(nC,:,nT)=smooth(data(nC,:,nT),smth);
                data(nC,:,nT)=data(nC,:,nT)-mean(data(nC,find(D.time<0),nT),2);
            end
        end
    end
    % LRP formula
    % LRP= 0.5 [meanleft hand response(C4V–C3V) + meanright hand response(C3V–C4V)]
    if avE
    LRP = 0.5 * (...
        mean(mean(data(C4_idx,:,leftResp_idx)-data(C3_idx,:,leftResp_idx),1),3) + ...
        mean(mean(data(C3_idx,:,rightResp_idx)-data(C4_idx,:,rightResp_idx),1),3));
    else
        LRP = 0.5 * squeeze(...
        mean(data(C4_idx,:,leftResp_idx)-data(C3_idx,:,leftResp_idx),3) + ...
        mean(data(C3_idx,:,rightResp_idx)-data(C4_idx,:,rightResp_idx),3));
    end
    
    if plotFlag
        figure;
        subplot 311
        plot(D.time,squeeze(mean(data(C4_idx,:,leftResp_idx),3)),'r-')
        hold on;
        plot(D.time,squeeze(mean(data(C3_idx,:,leftResp_idx),3)),'r--')
        plot(D.time,squeeze(mean(data(C3_idx,:,rightResp_idx),3)),'b-')
        plot(D.time,squeeze(mean(data(C4_idx,:,rightResp_idx),3)),'b--')
        legend({'rEl-lRe','lEl-lRe','lEl-rRe','rEl-rRe'})
        
        subplot 312
        plot(D.time,squeeze(mean(data(C4_idx,:,leftResp_idx),3))-squeeze(mean(data(C3_idx,:,leftResp_idx),3)),'r-')
        hold on;
        plot(D.time,squeeze(mean(data(C3_idx,:,rightResp_idx),3))-squeeze(mean(data(C4_idx,:,rightResp_idx),3)),'b-')
        legend({'lRe: contra-ipsi','rRe: contra-ipsi'})
        
        subplot 313
        plot(D.time, LRP,'LineWidth',2,'Color','r')
        
    end
end