function [p,c]=get_TPM(x,val)

if nargin<2
    val=[]; % expected values in r
end

% prepare
if isempty(val)
    m = unique(x); % get all levels data
else
    m=val;
end
n = length(x); % get size data
y = zeros(length(m),1);
p = zeros(length(m),length(m));
for k=1:n-1
    y(find(m==x(k))) = y(find(m==x(k))) + 1;
    p(find(m==x(k)),find(m==x(k+1))) = p(find(m==x(k)),find(m==x(k+1))) + 1;
end
c=p;
p = bsxfun(@rdivide,p,y); p(isnan(p)) = 0;
