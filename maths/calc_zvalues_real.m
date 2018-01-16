function zval = calc_zvalues_real(data,dim,mu0)

if nargin<3
    mu0 = 0;
end

data_mean = mean(data,dim);
data_std = std(data,[],dim);
n = size(data,dim);

zval = (data_mean-mu0)./(data_std);