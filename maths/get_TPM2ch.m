function [p,c]=get_TPM2ch(x,y)


% prepare
    m = unique(x); % get all levels data
n = length(x); % get size data
y = zeros(length(m),1);
p = zeros(2*length(m),2*length(m));
for k=1:n-1
    y(find(m==x(k))) = y(find(m==x(k))) + 1;
    p(find(m==x(k)),find(m==x(k+1))) = p(find(m==x(k)),find(m==x(k+1))) + 1;
end
c=p;
p = bsxfun(@rdivide,p,y); p(isnan(p)) = 0;
