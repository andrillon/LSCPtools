function [snr,f,logS,S]=get_power_SNR(data,params)

% params.tapers=[3 1];
% params.Fs=D.fsample;
% params.fpass = [1 55];
[S,f]=mtspectrumc(data,params);
logS=log(S);
kernel=params.kernel; %[-.25 -.25 0 0 1 0 0 -.25 -.25];
for nt=1:size(logS,2)
    snr(nt,:)= conv(logS(:,nt)', kernel, 'same');
end

logS=logS';
S=S';