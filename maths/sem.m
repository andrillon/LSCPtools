function SEM=sem(X,dim)
if nargin<2
    dim=1;
end
if isvector(X)
    SEM=nanstd(X)/sqrt(sum(~isnan(X))-1);
else
    SEM=nanstd(X,[],dim)/sqrt(min(sum(~isnan(X),dim))-1);
end