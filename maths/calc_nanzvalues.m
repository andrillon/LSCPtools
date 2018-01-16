function zval = calc_nanzvalues(data,dim,mu0)

if nargin<3
    mu0 = 0;
end

data_mean = nanmean(data,dim);
data_std = nanstd(data,[],dim);
n = size(data,dim);

zval = (data_mean-mu0)./(data_std/sqrt(n));

% error('Z-scores not calculated correctly - should be: (data-mu)/std')

