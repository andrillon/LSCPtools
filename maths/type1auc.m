%% 2015-12-17 Julian - Type-I AUC calculator
% Objective Performance modified for -4 to 4 confidence judgments
% Nonparametric estimate of task sensitivity
%
% Input signal truth ('signal') = [-1 1] or Left & Right
% Input choices ('decision') = [-4 : 4] or 1:4 Left & 1:4 Right
%
% Output is a pure frequency table and cumulative frequencies
% Also Type-I AUC calculated using AreaUnderROC
%
% Can also output d-prime and criterion

function [Frequencies, Cumulative_Frequencies, Type_One_AUC,dprime,criterion] = type1auc(signal,decision)

Frq = zeros(2,8);

for row = 1:length(signal)
    for confidence = 1:4 % Right responses
        if decision(row) == confidence && signal(row) > 0
            Frq(2,confidence+4) = (Frq(2,confidence+4)+1); % correct_rejection
        elseif decision(row) == confidence && signal(row) < 0
            Frq(1,(confidence+4)) = (Frq(1,(confidence+4)))+1; % miss
        end
    end
    for confidence = -4:-1 % Left responses
        if decision(row) == confidence && signal(row) < 0
            Frq(1,(confidence+5)) = (Frq(1,(confidence+5)))+1; % hit
        elseif decision(row) == confidence && signal(row) > 0
            Frq(2,(confidence+5)) = (Frq(2,(confidence+5)))+1; % false alarm
        end
    end
    
end

% Save Frequencies table
Frequencies = Frq;

% Convert to proportions
for column = 1:length(Frq)
    Frq_Pro(1,column) = (Frq(1,column)/sum(Frq(1,:)));
    Frq_Pro(2,column) = (Frq(2,column)/sum(Frq(2,:)));
end

% Cumulative proportions
for column = 1:length(Frq_Pro)
    Frq_Cum(1,column) = sum(Frq_Pro(1,(1:column)));
    Frq_Cum(2,column) = sum(Frq_Pro(2,(1:column)));
end

% Save Cumulative Frequencies table
Cumulative_Frequencies = Frq_Cum;

% Alternative method using trapz seems to give same result as AUC function
% Check the plot:
% plot([0 Frq_Cum(2,:)],[0 Frq_Cum(1,:)],'g:')
% trapz([0 Frq_Cum(2,:)],[0 Frq_Cum(1,:)])
% [X,Y,T,AUC]=perfcurve(signal,decision,1)

% Calculate area under ROC curve and save as Type_One_AUC
Type_One_AUC = AreaUnderROC([Frq_Cum(1,:);Frq_Cum(2,:)]');

if sum(Frq(1,5:8))==0 && sum(Frq(2,1:4))>0  && Type_One_AUC==1
    % Type 1 AUC == 1 but there are a >0 number of False Alarms
%     disp(['Perfect classification reported as a product of class imbalance'...
%         ', AUC recorded as NaN. Consider subject removal or use of d-prime.'])
    Type_One_AUC = NaN;
elseif sum(Frq(1,5:8))==0 && sum(Frq(2,1:4))>0
    % See for details: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4356897/
%     disp('No False Positives reported, AUC might be biased in this instance.')
end

%% Dprime & CRITERION: ASSUMES -1 IS HIT
% Includes loglinear correction approach from Hautus,1995
% Add 0.5 to hits & FAs, add 1 to Signal-Presence & Signal-Absence

% HITS / HITS + MISSES
hitrate = (sum(Frequencies(1,1:4))+.5)/(sum(Frequencies(1,:))+1);

% FALSE ALARMS / FALSE ALARMS + CORRECT REJECTIONS
falserate = (sum(Frequencies(2,1:4))+.5)/(sum(Frequencies(2,:))+1);

[dprime,criterion] = dprime_simple(hitrate,falserate);

end

%% INCLUDE AreaUnderROC FUNCTION TO MAKE PORTABLE

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

function [dp,c] = dprime_simple(h,fA)
% DPRIME_SIMPLE d' calculation given hit/false alarm rates
%   dp = dprime_simple(h,fA) returns the d' value dp 
%   for the hit rate h and false alarm rate fA
%   Karin Cox 8/31/14
%   updates:
%   8/31/14 added criterion output c, input arg checks   
%
%   inputs: 
%   h = hit rate, as float (0 < h < 1) 
%   fA = false alarm rate, as float (0 < fA < 1)
%
%   outputs: 
%   dp = d'
%   c = criterion c (negative values --> bias towards yes responses)
%
%   Example:
%   [dp,c] = dprime_simple(0.9,0.05)   
%   dp =
%     2.9264
%   c = 
%     0.1817
%
%   formulas: Stanislaw, H., & Todorov, N. (1999). Calculation of signal 
%   detection theory measures. Behavior research methods, instruments, & 
%   computers, 31(1), 137-149.


% check n args
narginchk(2,2);

% check for values out of bounds, also issue errors or warnings if =1 or 0
if or(or(h>1,h<0),or(fA>1,fA<0))
    error('input arguments must fall in the 0 to 1 range')
% standard d' formula returns NaN or Inf results for 0 or 1 inputs,
% corrections may be required (see article above for suggestions)
elseif or(or(h==1,h==0),or(fA==1,fA==0))
    warning('This function will not return finite values when h or fA = 0 or 1.')
end


% d prime = z(h)-z(fA)
dp = norminv(h)-norminv(fA);

% c = -0.5*[z(h)+z(fA)]
c = -0.5*(norminv(h)+ norminv(fA));

end
