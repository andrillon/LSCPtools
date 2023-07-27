function [real_out, perm_out]=lm_perm(table,predictor,formula,totperm,flagbinomial)
% table=GO_table;
%
% formula='RT~1+pred+Task+(1|SubID)';
% permvars={'Task','SubID'};
% totperm=100;
if nargin<5
    flagbinomial=0;
end
% run real model
eval(sprintf('table.pred=table.%s;',predictor));
if ~flagbinomial
    model= fitlme(table,formula);
else
    model= fitglme(table,formula,'Distribution','Binomial');
end
real_out=[double(model.Coefficients(4,2)) double(model.Coefficients(4,4)) double(model.Coefficients(4,6))];

perm_out=nan(totperm,4);
fprintf('%4.0f/%4.0f\n',0,totperm)
for np=1:totperm
    pred_perm=table.pred;
    pred_perm=pred_perm(randperm(length(pred_perm)));
    table2=table;
    table2.pred=pred_perm;
    if ~flagbinomial
        model= fitlme(table2,formula);
    else
        model= fitglme(table2,formula,'Distribution','Binomial');
    end
    perm_out(np,:)=[double(model.Coefficients(4,2)) double(model.Coefficients(4,4)) double(model.Coefficients(4,6)) np];
    fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',np,totperm)
end
fprintf('\n');