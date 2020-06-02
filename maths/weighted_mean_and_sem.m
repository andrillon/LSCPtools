function [M,S]=weighted_mean_and_sem(V,W)

W2=W;
V2=V;

W2(isnan(V) | W==0)=[];
V2(isnan(V) | W==0)=[];

M=nansum(V2.*W2)./sum(W2);
S=std(V2,W2)/sqrt(length(V2)-1);
   