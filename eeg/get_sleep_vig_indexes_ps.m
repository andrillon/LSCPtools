function [alpha_theta,vig_index,W_index,NREM_index,REM_index, powbyband]=get_sleep_vig_indexes_ps(pow,faxis,param)

% Set parameters
if ~isfield(param,'delta_band')
    param.delta_band=[0.1 2];
end
if ~isfield(param,'alpha_band')
    param.alpha_band=[8 13];
end
if ~isfield(param,'theta_band')
    param.theta_band=[4 7];
end
if ~isfield(param,'spindle_band')
    param.spindle_band=[11 16];
end
if ~isfield(param,'beta_band')
    param.beta_band=[20 40];
end

if isfield(param,'StopFreqs')
    StopFreqs=param.StopFreqs;
    StopFreqsIdx=[];
    for nfreq=1:length(StopFreqs)
        StopFreqsIdx=[StopFreqsIdx find(faxis>=StopFreqs(nfreq)-0.1 & faxis<=StopFreqs(nfreq)+0.1)];
    end
    pow(StopFreqsIdx)=[];
    faxis(StopFreqsIdx)=[];
end

alpha_pow=mean(pow(faxis>=param.alpha_band(1) & faxis<=param.alpha_band(2)));
theta_pow=mean(pow(faxis>=param.theta_band(1) & faxis<=param.theta_band(2)));
beta_pow=mean(pow(faxis>=param.beta_band(1) & faxis<=param.beta_band(2)));
spindle_pow=mean(pow(faxis>=param.spindle_band(1) & faxis<=param.spindle_band(2)));
delta_pow=mean(pow(faxis>=param.delta_band(1) & faxis<=param.delta_band(2)));
powbyband.alpha=alpha_pow;
powbyband.theta=theta_pow;
powbyband.delta=delta_pow;
powbyband.spindle=spindle_pow;
powbyband.beta=beta_pow;

% Extract alpha theta ratio
alpha_theta=alpha_pow./theta_pow;

% Vigilance index as defined in Kouider, Andrillon et al. Current Biology 2014
% VI = [delta power + theta power + spindle power] / [alpha power + high-beta power]
% With delta corresponding to 0.1 ? 4 Hz, theta to 4 ? 7 Hz, spindle frequency to 11 ? 16 Hz, alpha to 8 ? 13 Hz, and high-beta to 20 ? 40 Hz).
vig_index=(delta_pow+theta_pow+spindle_pow)./(alpha_pow+beta_pow);

% Wake Index
W_index=(alpha_pow+beta_pow)./(delta_pow+theta_pow+spindle_pow);

% NREM Index
NREM_index=(delta_pow+spindle_pow)./(alpha_pow+beta_pow);

% Wake Index
REM_index=(theta_pow+beta_pow)./(delta_pow+alpha_pow+spindle_pow);
