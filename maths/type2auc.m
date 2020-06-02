%% TYPE 2 AUC: Kunimoto & Perfcurve method
%
% Computes Type 2 AUC using MATLAB perfcurve function and fitglm.
% Input vectors of 'correct' and 'confidence' for default behaviour.
%
% correct = 1 or 0 ('correct' to 'incorrect')
% confidence = 1:4 (absolute 'low' to 'high' confidence respectively)
%
% Fits using Kunimoto technique if the maximum confidence judgement is <4
% (Perfcurve overestimates metacognition in this case) or Perfcurve_flag is
% specified. See Kaunitz et al (2016, PsychSci) for description of this.
%
% If only a single confidence level or correctness is recorded, by default
% a NaN is returned. Specify roc_flag to bypass NaN and compute
% metacognition using Fleming & Lau's (2014) type2roc.m technique (assuming
% it is added to path: see github.com/julian-matthews/stats-tools
%
% Julian Matthews 21/09/2017

function type2 = type2auc(correct,confidence,roc_flag,perfcurve_flag)

% Assume that perfcurve method is wanted unless flagged
if nargin < 3
    perfcurve_flag = 1;
    roc_flag = 0;
elseif nargin == 3
    perfcurve_flag = 1;
end

if isempty(confidence)
    % Empty input: can happen on iterative data processing
    disp('Empty input, coding as NaN');
    type2 = NaN;
else
    % Ensure absolute confidence values
    confidence = abs(confidence);
    
    %% CONTROL FOR MISSING CLASSIFIERS
    if isscalar(unique(correct)) && roc_flag == 0
        disp('All responses are either correct or incorrect, Type2 coded as NaN')
        type2 = NaN;
    elseif isscalar(unique(correct)) && roc_flag ~= 0
        disp('All responses are either correct or incorrect but roc_flagged. Type2roc employed with confid levels = 4')
        type2 = type2roc(correct,confidence,4); % Assumes 4-level confidence
    elseif isscalar(unique(confidence)) && roc_flag == 0
        disp('Only one confidence level specified, Type2 coded as NaN')
        type2 = NaN;
    elseif isscalar(unique(confidence)) && roc_flag ~= 0
        disp('Only one confidence level specified but roc_flagged. Type2roc employed with confid levels = 4')
        type2 = type2roc(correct,confidence,4); % Assumes 4-level confidence
    else
        
        clear type2 hitrate FArate
        
        if perfcurve_flag ~= 1 || max(confidence) < 4
            %% METHOD BY AreaUnderROC
            if perfcurve_flag == 1 && max(confidence) < 4
                disp('Confidence level 4 not used, employing Kunimoto technique')
            end
            
            for conf = 1:3
                
                % Compute number of metacognitive 'hits'
                rows = confidence>conf & correct==1;
                hits = sum(rows);
                
                % Compute number of metacognitive 'false alarms'
                rows = confidence>conf & correct==0;
                FAs = sum(rows);
                
                % Calculate hit rate and false alarm rate
                HR = hits / sum(correct==1); hitrate(conf) = HR; %#ok<*AGROW>
                FAR = FAs / sum(correct==0);  FArate(conf) = FAR;
                
            end
            
            type2 = AreaUnderROC([1 hitrate 0; 1 FArate 0]');
            
        else
            %% METHOD BY PERFCURVE
            
            resp = correct; % 1/0 for correct/incorrect
            pred = confidence; % Vector from 1:4
            
            mdl = fitglm(pred,resp,'Distribution','binomial',...
                'Link','logit','Intercept',false);
            scores = mdl.Fitted.Probability;
            
            [~,~,~,type2] = perfcurve(correct,scores,1);
            
        end
    end
end
end

%% Basic AreaUnderROC function for portability
function area = AreaUnderROC(ROC)
% area = AreaUnderROC(ROC)
%
% Compute the area under an ROC curve.
%
% xx/xx/xx  dhb  Wrote it.

[NFA,index] = sort(ROC(:,2));
NROC(:,1) = ROC(index,1);
NROC(:,2) = ROC(index,2);
[m,n] = size(NROC);

FROC = zeros(m+2,n);
FROC(1,:) = [0, 0];
FROC(m+2,:) = [1, 1];
FROC(2:m+1,:) = NROC;

area = 0;
totwidth = 0;
for i = 1:m+1
    meanhgt = (FROC(i,1)+FROC(i+1,1))/2;
    width = FROC(i+1,2)-FROC(i,2);
    area = area+meanhgt*width;
    totwidth = totwidth+width;
end

end