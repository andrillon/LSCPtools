function [d] = calc_dprime(Hit,FA)
%[d] = calc_dprime(Hit,FA)

% controle size
if size(Hit,1)>size(Hit,2)
    Hit=Hit';
end
if size(FA,1)>size(FA,2)
    FA=FA';
end

% take NaNs out
Hit(isnan(Hit))=[];
FA(isnan(FA))=[];

% control for Inf dprime
if sum(Hit==1)==0,
    Hit=[Hit 1];
elseif sum(Hit==0)==0,
    Hit=[Hit 0];
end
if sum(FA==1)==0,
    FA=[FA 1];
elseif sum(FA==0)==0,
    FA=[FA 0];
end

% compute dprime
d=(mean(Hit)-mean(FA))/...
    sqrt(0.5*(std(Hit)^2+std(FA)^2));
