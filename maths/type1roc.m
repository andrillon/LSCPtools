function [auroc1, cum_H1, cum_FA1] = type1roc(decision, stimulus, conf, Nratings)
% TA

i = Nratings+1;
for c = 1:Nratings
    H1(i-1) = length(find(conf == c & decision==1 & stimulus==1)) + 0.5;
    FA1(i-1) = length(find(conf == c & decision==1 & stimulus==0)) + 0.5;
    i = i-1;
end

H1save = H1; FA1save = FA1;

H1 = H1./sum(H1);
FA1 = FA1./sum(FA1);
cum_H1 = [0 cumsum(H1)];
cum_FA1 = [0 cumsum(FA1)];

i=1;
for c = 1:Nratings
    k(i) = (cum_H1(c+1) - cum_FA1(c))^2 - (cum_H1(c) - cum_FA1(c+1))^2;
    i = i+1;
end
auroc1 = 0.5 + 0.25*sum(k);