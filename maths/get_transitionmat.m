function p=get_transitionmat(x,states,corr)

m = length(states);
n = numel(x);
y = zeros(m,1);
p = zeros(m,m);
for k=1:n-1
    y(x(k)) = y(x(k)) + 1;
    p(x(k),x(k+1)) = p(x(k),x(k+1)) + 1;
end
p = bsxfun(@rdivide,p,y); p(isnan(p)) = 0;
if corr
for ns=1:length(states)
    if isempty(find(x==states(ns)))
        p(ns,:)=NaN;
        p(:,ns)=NaN;
    end
end
end